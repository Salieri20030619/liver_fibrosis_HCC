from pathlib import Path
import sys,json,pickle,hashlib,time,datetime,tarfile
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));N=P/'continuation_newnetwork_models_v2';I=N/'features';S=N/'search';V=P/'independent_validation_A_20261021_v1';O=P/'A_top5_external_20261021_v1';O.mkdir(exist_ok=False)
import pandas as pd,numpy as np,joblib
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def js(x,n):(O/n).write_text(json.dumps(x,indent=2,default=str))
t0=time.time();r=pd.read_csv(S/'all243_model_metrics.tsv',sep='\t');assert len(r)==243;r2=r.sort_values(['CV_mean_fold_AUC','model'],ascending=[False,True]);assert list(r.model)==list(r2.model);top=r.head(5).copy();top.insert(0,'training_rank',range(1,6));names=list(top.model);top[['training_rank','model','kind','CV_mean_fold_AUC','CV_fold_SD','CV_pooled_OOF_AUC']].to_csv(O/'frozen_training_top5.tsv',sep='\t',index=False)
with open(S/'meta_models.pkl','rb') as f:metas=pickle.load(f)
pairs={n:n.removeprefix('Stack_').split('__') for n in names};base_names=sorted(set(sum(pairs.values(),[])));assert len(base_names)==6;models={n:joblib.load(I/(n+'.joblib')) for n in base_names};assert all(n in metas for n in names)
files=[S/'all243_model_metrics.tsv',S/'frozen_search_rules.json',P/'cont_newnetwork_generated_combine.py',S/'meta_models.pkl',S/'all243_external_predictions.tsv',S/'all243_training_OOF.tsv',I/'train_features.tsv',I/'external_features.tsv',V/'new_features_original_scale.tsv']+[I/(n+'.joblib') for n in base_names];manifest={str(p):sha(p) for p in files};js(manifest,'source_hashes.json');js({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'models':names,'ranking':'Original saved outer10 meanAUC descending/modelname ascending','predictions':'Own saved base pipeline scaler/RFE and own meta; no fit/search/threshold adaptation','CI':'stratified percentile bootstrap2000seed20261021 sharedindices within cohort;differencevsrank1 pairedwithincohort;no crosscohortpairedCI','thresholds':'Only rank1 saved threshold confirmed; do not newly derive thresholds for othermodels','scope':'Exploratory external viewoftrainingselectedtop5;notreplacingbaseline'},'frozen_plan.json')
X=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0);E1=pd.read_csv(I/'external_features.tsv',sep='\t',index_col=0);E2=pd.read_csv(V/'new_features_original_scale.tsv',sep='\t',index_col=0);tm=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0).loc[X.index];m1=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0).loc[E1.index];m2=pd.read_csv(V/'new_sample_alignment.tsv',sep='\t',index_col=0).loc[E2.index];fold=pd.read_csv(I/'outer_cv_folds.tsv',sep='\t').set_index('sample').loc[X.index,'fold'];trainpred=pd.read_csv(S/'all243_training_OOF.tsv',sep='\t',index_col=0).loc[X.index,names];old=pd.read_csv(S/'all243_external_predictions.tsv',sep='\t',index_col=0).loc[E1.index,names]
for E in [E1,E2]:assert list(E.columns)==list(X.columns) and np.isfinite(E).all().all()
checks=[];features=[]
for n,m in models.items():
 assert m.n_features_in_==30 and list(m.classes_)==[0,1];assert list(m.named_steps)==['scale','RFE','model'];mask=m.named_steps['RFE'].support_;assert len(mask)==30 and mask.sum()==15
 for gene,used in zip(X.columns,mask):features.append(dict(base_model=n,gene=gene,used=bool(used),GSE162694_available=True,GSE174478_available=True))
 checks.append(dict(base_model=n,model_path=str(I/(n+'.joblib')),sha256=sha(I/(n+'.joblib')),input_features=30,selected_features=int(mask.sum()),own_scaler=True,own_RFE=True))
pd.DataFrame(features).to_csv(O/'own_base_feature_masks.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(O/'model_cache_checks.tsv',sep='\t',index=False)
for z in top.itertuples():
 aucs=[roc_auc_score(tm.y[fold==k],trainpred.loc[fold==k,z.model]) for k in range(10)];assert abs(np.mean(aucs)-z.CV_mean_fold_AUC)<1e-12 and abs(roc_auc_score(tm.y,trainpred[z.model])-z.CV_pooled_OOF_AUC)<1e-12
metrics=[];diff=[];replay=[]
for cohort,E,m in [('GSE162694',E1,m1),('GSE174478',E2,m2)]:
 bp=pd.DataFrame({n:models[n].predict_proba(E.to_numpy())[:,1] for n in base_names},index=E.index);pred=pd.DataFrame(index=E.index)
 for n in names:
  assert list(metas[n].feature_names_in_)==pairs[n];pred[n]=metas[n].predict_proba(bp[pairs[n]])[:,1]
  if cohort=='GSE162694':err=float(np.max(abs(pred[n]-old[n])));assert err<1e-12
  else:err=float(np.max(abs(pred[n]-metas[n].predict_proba(bp[pairs[n]])[:,1])))
  replay.append(dict(cohort=cohort,model=n,max_prediction_error=err))
 pred.to_csv(O/(cohort+'_predictions_LOCAL_ONLY.tsv'),sep='\t');y=m.y.to_numpy();rng=np.random.default_rng(20261021);pos=np.flatnonzero(y==1);neg=np.flatnonzero(y==0);ix=np.array([np.r_[rng.choice(pos,len(pos)),rng.choice(neg,len(neg))] for _ in range(2000)]);assert np.array_equal(ix,np.load(V/(cohort+'_bootstrap_indices.npy')));np.save(O/(cohort+'_bootstrap_indices.npy'),ix);boot=np.array([[roc_auc_score(y[ii],pred[n].to_numpy()[ii]) for n in names] for ii in ix]);np.save(O/(cohort+'_bootstrap_AUC.npy'),boot)
 for j,n in enumerate(names):
  lo,hi=np.quantile(boot[:,j],[.025,.975]);auc=roc_auc_score(y,pred[n]);metrics.append(dict(training_rank=j+1,model=n,cohort=cohort,cases=len(pos),controls=len(neg),AUC=auc,CI_low=lo,CI_high=hi,includes_0_5=bool(lo<=.5<=hi),CI_position='above_0.5' if lo>.5 else ('below_0.5' if hi<.5 else 'includes_0.5'),AP=average_precision_score(y,pred[n]),Brier=brier_score_loss(y,pred[n]),bootstrap=2000,seed=20261021))
  if j:
   dl,dh=np.quantile(boot[:,j]-boot[:,0],[.025,.975]);diff.append(dict(training_rank=j+1,model=n,reference=names[0],cohort=cohort,AUC_difference=auc-roc_auc_score(y,pred[names[0]]),CI_low=dl,CI_high=dh,includes_zero=bool(dl<=0<=dh)))
pd.DataFrame(metrics).sort_values(['training_rank','cohort']).to_csv(O/'top5_external_AUC_CI.tsv',sep='\t',index=False);pd.DataFrame(diff).to_csv(O/'paired_difference_vs_rank1.tsv',sep='\t',index=False);pd.DataFrame(replay).to_csv(O/'prediction_replay.tsv',sep='\t',index=False);assert all(sha(Path(p))==h for p,h in manifest.items());js({'status':'completed','fit_calls':0,'base_models_loaded':len(models),'meta_models_used':5,'own_preprocessing':True,'coverage':'30/30allinputand15/15eachbase;noimputation','hashes_unchanged':True,'training_ranking_replayed':True,'max_prediction_error':max(z['max_prediction_error'] for z in replay),'bootstrap_indices_identical_to_prior_validation':True,'seconds':time.time()-t0},'validation.json');print(pd.DataFrame(metrics).to_string(index=False),flush=True)
(O/'top5_worker_20261021.py').write_bytes(Path(__file__).read_bytes());files=[p for p in O.iterdir() if p.is_file()];pd.DataFrame([dict(path=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in files]).to_csv(O/'artifact_manifest.tsv',sep='\t',index=False);stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');arc=P/('top5_delivery_'+stamp+'.tar.gz')
with tarfile.open(arc,'w:gz') as t:
 for p in files+[O/'artifact_manifest.tsv']:t.add(p,arcname='evidence/'+p.name)
print(json.dumps({'archive':str(arc),'bytes':arc.stat().st_size,'sha256':sha(arc),'utc':stamp}),flush=True)
