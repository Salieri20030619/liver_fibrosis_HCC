from pathlib import Path
import os,sys,json,hashlib,time,traceback,itertools,pickle
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));O=P/'candidate_ablation_20261021_v1';N=P/'continuation_newnetwork_models_v2';I=N/'features';S=N/'search'
import numpy as np,pandas as pd,joblib
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def js(d,p):p.write_text(json.dumps(d,indent=2,default=str))
def dump(m,p):tmp=p.with_suffix('.tmp');joblib.dump(m,tmp);os.replace(tmp,p)
tm=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0);baseX=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0).loc[tm.index];fold=pd.read_csv(I/'outer_cv_folds.tsv',sep='\t').set_index('sample').loc[tm.index,'fold'];y=tm.y
genes50=pd.read_csv(O/'RSS_union50.tsv',sep='\t').gene.tolist();assert set(baseX.columns)<=set(genes50)
alltrain=pd.read_csv(P/'F6_S4_S6_F8_gene_validation/GSE135251_log2CPM_full.tsv.gz',sep='\t',index_col=0);wide=alltrain.loc[genes50,tm.index].T;assert np.isfinite(wide.to_numpy()).all()
arms={'A_RFE30':baseX,'B_all30':baseX,'C_RSS50':wide}
plan={'created_UTC':pd.Timestamp.now(tz='UTC').isoformat(),'status':'frozen_before_ablation_fits','arms':{'A_RFE30':'Historical30 features, train-fold RFE retain15; replay fixed baseline only','B_all30':'Same30, no RFE','C_RSS50':'All existing six-state top10RSS union50; remove only CSI-module/highest-state intersection; no RFE'},'C_choice':'Minimal proposed extension; alternative complete329 regulons disclosed but not run. No new topN or p threshold. 50 original TFs all mapped in training.','comparison':'A vs B isolates RFE; B vs C isolates candidate pool at identical no-RFE pipeline','model':'RF200 trees depth4 minleaf3 max_features=sqrt random_state3; SVM RBF C1 gamma=scale probability=True random_state3; logistic meta C1 maxiter2000','tuning_budget':'one prespecified configuration per arm, no hyperparameter or243 search','splits':'Exact historical outer10 and inner5 index files; inner5 solely cross-fitting meta inputs','within_training':'Scaler and A RFE fit inside every base training partition; no label-driven screening on alltrain; meta fit to innerOOF only','external':'Previously seen112 only after training ranking frozen; technical OVOL2 absence filled fulltraining mean, no gene removal/no external label optimization','selection':'rank by mean outer-foldAUROC; tie pooledBrier then armname; comparison exploratory, baseline not replaced','cache':'Each pair of base fits and predictions plus input fingerprints and memberships persisted; reload only if fingerprint matches; atomic joblib files. All newly fitted meta/base models retained.','no_new_SHAP':True}
planpath=O/'frozen_comparison_plan.json'
if planpath.exists():
 old=json.loads(planpath.read_text());assert old['arms']==plan['arms'] and old['model']==plan['model']
else:js(plan,planpath)
source=[I/'train_features.tsv',I/'train_design.tsv',I/'external_features.tsv',I/'external_design.tsv',I/'outer_cv_folds.tsv',I/'RF_4.joblib',I/'SVM_1.0.joblib',S/'meta_models.pkl',P/'continuation_SCENIC_candidates_v1/top10_RSS_per_state.tsv',P/'continuation_SCENIC_candidates_v1/F5G_candidates.tsv',P/'SHAP_repair_243_20261021_v1/validation.json']
manifest=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in source]
if (O/'baseline_source_manifest.tsv').exists():
 prior=pd.read_csv(O/'baseline_source_manifest.tsv',sep='\t');assert all(sha(Path(r.path))==r.sha256 for r in prior.itertuples())
else:pd.DataFrame(manifest).to_csv(O/'baseline_source_manifest.tsv',sep='\t',index=False)
for name,x in arms.items():pd.DataFrame({'order':range(1,len(x.columns)+1),'gene':x.columns}).to_csv(O/(name+'_frozen_features.tsv'),sep='\t',index=False)
fitlog=[];selection=[];partitionaudit=[];failures=[];preds={};metafits={};foldmetrics=[];start=time.time()
def basefit(arm,context,part,x,trainids,testids):
 q=O/'cache'/arm/context/part;q.mkdir(parents=True,exist_ok=True)
 sig=hashlib.sha256(x.loc[trainids].to_csv().encode()+y.loc[trainids].to_csv().encode()+x.loc[testids].to_csv().encode()+sha(Path(__file__)).encode()+json.dumps(plan['model']).encode()+arm.encode()).hexdigest()
 cached=(q/'complete.json').exists()
 if cached:assert json.loads((q/'complete.json').read_text())['signature']==sig
 result={};models={}
 for name in ['RF_4','SVM_1.0']:
  t=time.time();path=q/(name+'.joblib')
  receipt=q/(name+'_complete.json')
  if receipt.exists():
   assert json.loads(receipt.read_text())['signature']==sig
   m=joblib.load(path);action='loaded_cached'
  else:
   steps=[('scale',StandardScaler())]
   if arm=='A_RFE30':steps.append(('RFE',RFE(LogisticRegression(C=1,solver='liblinear',max_iter=2000),n_features_to_select=15,step=1)))
   est=RandomForestClassifier(n_estimators=200,max_depth=4,min_samples_leaf=3,n_jobs=4,random_state=3) if name=='RF_4' else SVC(C=1,probability=True,random_state=3)
   m=Pipeline(steps+[('model',est)]);m.fit(x.loc[trainids].to_numpy(),y.loc[trainids]);dump(m,path);js({'signature':sig},receipt);action='fit_new'
  assert list(m.classes_)==[0,1]
  models[name]=m;result[name]=m.predict_proba(x.loc[testids].to_numpy())[:,1] if len(testids) else np.array([])
  features=x.columns[m.named_steps['RFE'].support_].tolist() if arm=='A_RFE30' else x.columns.tolist()
  selection.extend(dict(arm=arm,context=context,partition=part,base_model=name,gene=g,selected=g in features) for g in x.columns)
  fitlog.append(dict(arm=arm,context=context,partition=part,model=name,action=action,seconds=time.time()-t,train_n=len(trainids),test_n=len(testids),features=len(features),path=str(path)))
 pd.DataFrame({'sample':trainids,'role':'train'}).to_csv(q/'training_members.tsv',sep='\t',index=False);pd.DataFrame({'sample':testids,'role':'test'}).to_csv(q/'heldout_members.tsv',sep='\t',index=False)
 out=pd.DataFrame(result,index=testids);out.to_csv(q/'heldout_base_predictions.tsv',sep='\t');js({'signature':sig,'status':'completed'},q/'complete.json');return out,models
for arm,x in arms.items():
 try:
  oof=pd.Series(index=x.index,dtype=float)
  for context in [str(k) for k in range(10)]+['full']:
   ids=x.index if context=='full' else x.index[fold!=int(context)];testids=x.index[:0] if context=='full' else x.index[fold==int(context)]
   inner=pd.read_csv(S/('inner_folds_'+context+'.tsv'),sep='\t').set_index('sample').loc[ids,'inner_fold'];assert set(inner.index)==set(ids)
   innerpred=pd.DataFrame(index=ids,columns=['RF_4','SVM_1.0'],dtype=float)
   for k in range(5):
    tr=ids[inner!=k];te=ids[inner==k];assert not set(tr)&set(te) and not set(testids)&set(tr) and not set(testids)&set(te)
    partitionaudit.append(dict(arm=arm,context=context,inner_fold=k,train_n=len(tr),heldout_n=len(te),outer_holdout_overlap=0))
    z,_=basefit(arm,context,'inner_'+str(k),x,tr,te);innerpred.loc[te]=z
   assert np.isfinite(innerpred.to_numpy()).all()
   q=O/'cache'/arm/context;innerpred.to_csv(q/'innerOOF_meta_training.tsv',sep='\t')
   mp=q/'meta.joblib'
   msig=hashlib.sha256(innerpred.to_csv().encode()+y.loc[ids].to_csv().encode()).hexdigest();mr=q/'meta_complete.json'
   if mr.exists():
    assert json.loads(mr.read_text())['signature']==msig
    meta=joblib.load(mp)
   else:
    meta=LogisticRegression(max_iter=2000).fit(innerpred,y.loc[ids]);dump(meta,mp);js({'signature':msig},mr)
   z,models=basefit(arm,context,'base_full',x,ids,testids)
   if context!='full':
    pp=meta.predict_proba(z)[:,1];oof.loc[testids]=pp;foldmetrics.append(dict(arm=arm,fold=int(context),n=len(testids),AUROC=roc_auc_score(y.loc[testids],pp),AP=average_precision_score(y.loc[testids],pp),Brier=brier_score_loss(y.loc[testids],pp)))
   else:metafits[arm]=(meta,models)
   pd.DataFrame(fitlog).to_csv(O/'fit_or_load_log.tsv',sep='\t',index=False);pd.DataFrame(selection).to_csv(O/'fold_feature_selections.tsv',sep='\t',index=False);pd.DataFrame(partitionaudit).to_csv(O/'partition_leakage_audit.tsv',sep='\t',index=False)
   js({'arm':arm,'context_completed':context,'elapsed_seconds':time.time()-start},O/'progress.json');print(arm,context,'completed',round(time.time()-start),flush=True)
  assert oof.notna().all();preds[arm]=oof;oof.to_csv(O/(arm+'_training_OOF.tsv'),sep='\t',header=['probability'])
 except Exception:
  failures.append({'arm':arm,'context':context,'error':traceback.format_exc()});js(failures,O/'failures.json');print('FAILED',arm,flush=True)
pd.DataFrame(foldmetrics).to_csv(O/'outer_fold_metrics.tsv',sep='\t',index=False)
js(failures,O/'failures.json')
assert len(preds)==3,'Incomplete arms explicitly retained; no partial mean comparison'
old=pd.read_csv(S/'all243_training_OOF.tsv',sep='\t',index_col=0).loc[tm.index,'Stack_RF_4__SVM_1.0'];aerr=float(np.max(abs(preds['A_RFE30']-old)));assert aerr<1e-9,aerr
rows=[]
for arm,pp in preds.items():
 fs=pd.DataFrame(foldmetrics).query('arm==@arm');rows.append(dict(arm=arm,mean_fold_AUROC=fs.AUROC.mean(),SD_fold_AUROC=fs.AUROC.std(),pooled_OOF_AUROC=roc_auc_score(y,pp),OOF_AP=average_precision_score(y,pp),OOF_Brier=brier_score_loss(y,pp)))
r=pd.DataFrame(rows).sort_values(['mean_fold_AUROC','OOF_Brier','arm'],ascending=[False,True,True]);r.to_csv(O/'training_internal_comparison.tsv',sep='\t',index=False)
js({'ranking':r.arm.tolist(),'selected_by_training_only':r.iloc[0].arm,'main_baseline_replaced':False,'baseline_OOF_replay_max_error':aerr,'external_labels_not_used':True},O/'TRAINING_COMPARISON_FROZEN.json')
# Only now read external labels and generate exploratory predictions.
em=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0);allE=pd.read_csv(P/'F6_S4_S6_F8_gene_validation/GSE162694_log2CPM_full.tsv.gz',sep='\t',index_col=0);external={};missing=[]
for arm,x in arms.items():
 e=allE.reindex(x.columns).loc[:,em.index].T
 for g in e:
  if e[g].isna().any():missing.append(dict(arm=arm,gene=g,n_missing=int(e[g].isna().sum()),replacement_training_mean=float(x[g].mean())));e[g]=e[g].fillna(x[g].mean())
 meta,models=metafits[arm];z=pd.DataFrame({n:m.predict_proba(e.to_numpy())[:,1] for n,m in models.items()},index=e.index);pp=meta.predict_proba(z)[:,1];external[arm]=pp
 extrow={'arm':arm,'AUROC':roc_auc_score(em.y,pp),'AP':average_precision_score(em.y,pp),'Brier':brier_score_loss(em.y,pp),'n':len(pp)}
 rows=[] if arm=='A_RFE30' else rows;rows.append(extrow)
pd.DataFrame(rows).to_csv(O/'external_exploratory_comparison.tsv',sep='\t',index=False);pd.DataFrame(external,index=em.index).to_csv(O/'external_predictions.tsv',sep='\t');pd.DataFrame(missing,columns=['arm','gene','n_missing','replacement_training_mean']).to_csv(O/'external_missing_feature_handling.tsv',sep='\t',index=False)
oldext=pd.read_csv(S/'all243_external_predictions.tsv',sep='\t',index_col=0).loc[em.index,'Stack_RF_4__SVM_1.0'];eerr=float(np.max(abs(external['A_RFE30']-oldext)));assert eerr<1e-9
# Paired bootstrap is conditional on the fixed OOF/external predictions, not a retraining CI.
boots=[];rng=np.random.default_rng(20261021)
for cohort,truth,pr in [('training_OOF',y.to_numpy(),{k:p.to_numpy() for k,p in preds.items()}),('external_exploratory',em.y.to_numpy(),external)]:
 for left,right in [('B_all30','A_RFE30'),('C_RSS50','B_all30'),('C_RSS50','A_RFE30')]:
  ds=[];pos=np.flatnonzero(truth==1);neg=np.flatnonzero(truth==0)
  for _ in range(2000):
   ii=np.r_[rng.choice(pos,len(pos)),rng.choice(neg,len(neg))];ds.append(roc_auc_score(truth[ii],pr[left][ii])-roc_auc_score(truth[ii],pr[right][ii]))
  boots.append(dict(cohort=cohort,comparison=left+' minus '+right,AUROC_difference=roc_auc_score(truth,pr[left])-roc_auc_score(truth,pr[right]),CI_low=np.quantile(ds,.025),CI_high=np.quantile(ds,.975),scope='conditional paired stratified bootstrap2000; no model refit; external not selection'))
pd.DataFrame(boots).to_csv(O/'paired_AUROC_comparisons.tsv',sep='\t',index=False)
sel=pd.DataFrame(selection);sf=sel[(sel.partition=='base_full')&(sel.context!='full')&(sel.base_model=='RF_4')];sf.groupby(['arm','gene']).selected.agg(['sum','mean','count']).reset_index().to_csv(O/'outer_feature_selection_stability.tsv',sep='\t',index=False)
assert all(sha(Path(row['path']))==row['sha256'] for row in manifest)
js({'status':'completed','arms':3,'outer_folds':10,'inner_folds':5,'baseline_OOF_error':aerr,'baseline_external_error':eerr,'source_hashes_unchanged':True,'failed_arms':len(failures),'fitted_base_models':sum(r['action']=='fit_new' for r in fitlog),'loaded_base_models':sum(r['action']=='loaded_cached' for r in fitlog),'meta_models':33,'partition_checks':len(partitionaudit),'training_only_choice':r.iloc[0].arm if isinstance(r,pd.DataFrame) else None,'elapsed_seconds':time.time()-start,'C_external_missing':['OVOL2'],'baseline_preserved':True},O/'validation.json')
print((O/'validation.json').read_text(),flush=True)
