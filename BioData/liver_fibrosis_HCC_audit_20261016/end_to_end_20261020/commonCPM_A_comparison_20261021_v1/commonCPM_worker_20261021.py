from pathlib import Path
import sys,os,json,time,hashlib,traceback,itertools
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));O=P/'commonCPM_A_comparison_20261021_v1';O.mkdir(exist_ok=True);B=P/'candidate_ablation_20261021_v1';I=P/'continuation_newnetwork_models_v2/features';S=P/'continuation_newnetwork_models_v2/search'
import numpy as np,pandas as pd,joblib
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest,f_classif,RFE
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def js(x,p):p.write_text(json.dumps(x,indent=2,default=str))
def dump(x,p):q=p.with_suffix('.tmp');joblib.dump(x,q);os.replace(q,p)
X=pd.read_csv(O/'train_features_commonCPM.tsv',sep='\t',index_col=0);genes=X.columns.tolist();tm=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0).loc[X.index];y=tm.y;fold=pd.read_csv(I/'outer_cv_folds.tsv',sep='\t').set_index('sample').loc[X.index,'fold'];assert X.shape==(206,30) and np.isfinite(X).all().all()
plan={'change':'F vs original A changes only CPM denominator from cohort full gene set to identity-only common30318; original30 candidates, fold-refit scaler/RFE15 and RF4/SVM1 stacking unchanged. No D/E rules.','common_gene_SHA256':sha(O/'common_Ensembl_30318.tsv'),'availability':'Known external platform gene availability used only as identity intersection; exploratory preprocessing; no external expression, labels or AUC choose denominator genes','algorithm':'RF200 depth4 minleaf3 n_jobs4 seed3; SVC C1 RBF gamma=scale probability=True seed3; RFE Logistic C1 liblinear maxiter2000 step1 retain15; meta Logistic C1 maxiter2000','validation':'Same outer10 inner5 for stacking, no tuning or243 search; A reused verified cache','CI':'stratified percentile bootstrap2000 seed20261021 shared resamples saved per cohort, fixed predictions','script_sha':sha(Path(__file__))}
if (O/'frozen_plan.json').exists():assert json.loads((O/'frozen_plan.json').read_text())==plan
else:js(plan,O/'frozen_plan.json')
pd.DataFrame({'gene':genes}).to_csv(O/'original_candidate30.tsv',sep='\t',index=False)
source=[O/'common_Ensembl_30318.tsv',O/'train_features_commonCPM.tsv',O/'external_features_commonCPM.tsv',I/'RF_4.joblib',I/'SVM_1.0.joblib',S/'meta_models.pkl',I/'outer_cv_folds.tsv',P/'SHAP_repair_243_20261021_v1/validation.json',B/'A_RFE30_training_OOF.tsv',B/'external_predictions.tsv'];manifest={str(p):sha(p) for p in source};js(manifest,O/'source_manifest.json')
logs=[];selections=[];scores=[];audit=[];start=time.time()
def fitpart(context,part,tr,te):
 q=O/'cache'/context/part;q.mkdir(parents=True,exist_ok=True);sig=hashlib.sha256(X.loc[tr].to_csv().encode()+y.loc[tr].to_csv().encode()+X.loc[te].to_csv().encode()+plan['script_sha'].encode()).hexdigest();out={};models={}
 for name in ['RF_4','SVM_1.0']:
  path=q/(name+'.joblib');receipt=q/(name+'_complete.json');t=time.time()
  if receipt.exists():
   assert json.loads(receipt.read_text())['signature']==sig;m=joblib.load(path);action='loaded'
  else:
   est=RandomForestClassifier(n_estimators=200,max_depth=4,min_samples_leaf=3,n_jobs=4,random_state=3) if name=='RF_4' else SVC(C=1,probability=True,random_state=3)
   m=Pipeline([('scale',StandardScaler()),('RFE',RFE(LogisticRegression(C=1,solver='liblinear',max_iter=2000),n_features_to_select=15,step=1)),('model',est)])
   m.fit(X.loc[tr].to_numpy(),y.loc[tr]);dump(m,path);js({'signature':sig},receipt);action='fit_new'
  selected=np.array(genes);final=selected[m.named_steps['RFE'].support_];assert len(selected)==30 and len(final)==15
  for g in genes:selections.append(dict(context=context,partition=part,model=name,gene=g,candidate30=g in selected,final15=g in final))
  assert np.max(abs(m.named_steps['scale'].mean_-X.loc[tr,selected].mean().to_numpy()))<1e-10
  out[name]=m.predict_proba(X.loc[te].to_numpy())[:,1] if len(te) else [];models[name]=m;logs.append(dict(context=context,partition=part,model=name,action=action,seconds=time.time()-t,train_n=len(tr),test_n=len(te),path=str(path)))
 pd.DataFrame({'sample':tr}).to_csv(q/'training_members.tsv',sep='\t',index=False);pd.DataFrame({'sample':te}).to_csv(q/'heldout_members.tsv',sep='\t',index=False);z=pd.DataFrame(out,index=te);z.to_csv(q/'heldout_base_predictions.tsv',sep='\t');return z,models
try:
 oof=pd.Series(index=X.index,dtype=float)
 for context in [str(k) for k in range(10)]+['full']:
  ids=X.index if context=='full' else X.index[fold!=int(context)];te=X.index[:0] if context=='full' else X.index[fold==int(context)];inner=pd.read_csv(S/('inner_folds_'+context+'.tsv'),sep='\t').set_index('sample').loc[ids,'inner_fold'];ip=pd.DataFrame(index=ids,columns=['RF_4','SVM_1.0'],dtype=float)
  for k in range(5):
   tr=ids[inner!=k];it=ids[inner==k];assert not set(tr)&set(it) and not set(te)&(set(tr)|set(it));ip.loc[it],_=fitpart(context,'inner_'+str(k),tr,it);audit.append(dict(context=context,inner=k,overlap=0))
  assert np.isfinite(ip).all().all();q=O/'cache'/context;ip.to_csv(q/'innerOOF_meta_training.tsv',sep='\t');mp=q/'meta.joblib';mr=q/'meta_complete.json';sig=hashlib.sha256(ip.to_csv().encode()+y.loc[ids].to_csv().encode()).hexdigest()
  if mr.exists():assert json.loads(mr.read_text())['signature']==sig;meta=joblib.load(mp)
  else:meta=LogisticRegression(max_iter=2000).fit(ip,y.loc[ids]);dump(meta,mp);js({'signature':sig},mr)
  z,models=fitpart(context,'base_full',ids,te)
  if context!='full':oof.loc[te]=meta.predict_proba(z)[:,1]
  pd.DataFrame(logs).to_csv(O/'fit_log.tsv',sep='\t',index=False);pd.DataFrame(selections).to_csv(O/'feature_selections.tsv',sep='\t',index=False);pd.DataFrame(scores).to_csv(O/'partition_F_scores.tsv',sep='\t',index=False);pd.DataFrame(audit).to_csv(O/'partition_audit.tsv',sep='\t',index=False);print(context,'completed',round(time.time()-start,2),flush=True)
 assert oof.notna().all();old=pd.read_csv(B/'A_RFE30_training_OOF.tsv',sep='\t',index_col=0).loc[X.index].iloc[:,0];pr=pd.DataFrame({'A_original':old,'F_commonCPM':oof});pr.to_csv(O/'training_OOF.tsv',sep='\t');js({'training_completed_before_external':True,'baseline_not_replaced':True},O/'TRAINING_FROZEN.json')
 em=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0);ext=pd.read_csv(O/'external_features_commonCPM.tsv',sep='\t',index_col=0).loc[em.index,genes];missing=[];assert np.isfinite(ext).all().all();z=pd.DataFrame({n:m.predict_proba(ext.to_numpy())[:,1] for n,m in models.items()},index=ext.index);ep=pd.DataFrame({'A_original':pd.read_csv(B/'external_predictions.tsv',sep='\t',index_col=0).loc[em.index,'A_RFE30'],'F_commonCPM':meta.predict_proba(z)[:,1]});ep.to_csv(O/'external_predictions.tsv',sep='\t')
 metrics=[];cis=[];diff=[];fm=[]
 for cohort,truth,pred in [('training_OOF',y.to_numpy(),pr),('external_exploratory',em.y.to_numpy(),ep)]:
  rng=np.random.default_rng(20261021);pos=np.flatnonzero(truth==1);neg=np.flatnonzero(truth==0);bb=[];indices=[]
  for j in range(2000):
   ix=np.r_[rng.choice(pos,len(pos)),rng.choice(neg,len(neg))];indices.append(ix);bb.append([roc_auc_score(truth[ix],pred[a].to_numpy()[ix]) for a in pred])
  bb=np.array(bb);np.save(O/(cohort+'_bootstrap_indices.npy'),np.array(indices))
  for j,a in enumerate(pred):
   auc=roc_auc_score(truth,pred[a]);lo,hi=np.quantile(bb[:,j],[.025,.975]);metrics.append(dict(cohort=cohort,arm=a,AUC=auc,CI_low=lo,CI_high=hi,includes_0_5=bool(lo<=.5<=hi),AP=average_precision_score(truth,pred[a]),Brier=brier_score_loss(truth,pred[a]),cases=len(pos),controls=len(neg),bootstrap=2000,seed=20261021))
  lo,hi=np.quantile(bb[:,1]-bb[:,0],[.025,.975]);diff.append(dict(cohort=cohort,comparison='F minus A',AUC_difference=roc_auc_score(truth,pred.iloc[:,1])-roc_auc_score(truth,pred.iloc[:,0]),CI_low=lo,CI_high=hi))
 for a in pr:
  for k in range(10):
   ii=fold==k;fm.append(dict(arm=a,fold=k,n=int(ii.sum()),AUC=roc_auc_score(y[ii],pr.loc[ii,a]),AP=average_precision_score(y[ii],pr.loc[ii,a]),Brier=brier_score_loss(y[ii],pr.loc[ii,a])))
 pd.DataFrame(metrics).to_csv(O/'AUC_CI_metrics.tsv',sep='\t',index=False);pd.DataFrame(diff).to_csv(O/'paired_differences.tsv',sep='\t',index=False);pd.DataFrame(fm).to_csv(O/'fold_metrics.tsv',sep='\t',index=False)
 sel=pd.DataFrame(selections);ss=sel[(sel.context!='full')&(sel.partition=='base_full')&(sel.model=='RF_4')];ss.groupby('gene')[['candidate30','final15']].sum().to_csv(O/'feature_stability.tsv',sep='\t');sel[(sel.context=='full')&(sel.partition=='base_full')&(sel.model=='RF_4')].to_csv(O/'final_selected_features.tsv',sep='\t',index=False)
 assert all(sha(Path(p))==h for p,h in manifest.items());js([],O/'failures.json');v={'status':'completed','seconds':time.time()-start,'new_base_fits':sum(r['action']=='fit_new' for r in logs),'cache_hits':sum(r['action']=='loaded' for r in logs),'meta_models':11,'baseline_reused_no_refit':True,'source_hashes_unchanged':True,'external_missing':missing,'outer_folds':10,'inner_folds':5,'partition_checks':len(audit)};js(v,O/'validation.json');print(json.dumps(v),flush=True)
except Exception:
 js([traceback.format_exc()],O/'failures.json');raise
