from pathlib import Path
import sys,os,time,json,pickle,subprocess
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));sys.path.insert(0,str(P/'shap_compatible'));D=P.parent/'analysis_stage1_20261017_075157UTC';N=P/'continuation_newnetwork_models_v2';I=N/'features'
for i in range(1440):
 if (N/'validation.json').exists():break
 time.sleep(10)
else:raise RuntimeError('Expandednetworkmodels pending4hours')
state=json.loads((N/'validation.json').read_text())
if state['status']=='reused_identical_feature_matrix':print('Reuse existing verifiedSHAP becauseavailablefeaturesunchanged');raise SystemExit(0)
assert state['status']=='models_and_metrics_complete_SHAP_pending',state
import numpy as np,pandas as pd,joblib,shap
O=N/'SHAP';O.mkdir(exist_ok=False);X=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0);E=pd.read_csv(I/'external_features.tsv',sep='\t',index_col=0);selected=json.loads((N/'search/validation.json').read_text())['selected']['model'];names=selected[len('Stack_'):].split('__') if selected.startswith('Stack_') else [selected];rnames=[n for n in names if n.startswith(('Stepglm','glmBoost','plsRglm'))];pymodels={n:joblib.load(I/(n+'.joblib')) for n in names if n not in rnames};meta=None
if selected.startswith('Stack_'):
 with open(N/'search/meta_models.pkl','rb') as f:meta=pickle.load(f)[selected]
proc=None;calls=0
if rnames:
 r=r'''.libPaths(c('LIB','OLD/irGSEA_dependencies_v1/R_project_library_v1','OLD/trajectory_dependencies_v1/R_project_library_v1',.libPaths()));library(mboost);library(plsRglm)
 names<-NAMES;models<-lapply(names,function(n)readRDS(paste0('MODELS/',n,'.rds')));names(models)<-names;con<-file('stdin',open='r');cat('READY\n');flush.console()
 repeat{line<-readLines(con,n=1);if(!length(line)||line=='STOP')break;x<-as.matrix(read.delim(line,check.names=FALSE));out<-list();for(n in names){m<-models[[n]];z<-as.data.frame(scale(x,m$preprocessing$mean,m$preprocessing$sd)[,m$preprocessing$features,drop=FALSE]);out[[n]]<-if(grepl('plsRglm',n))as.numeric(predict(m$fit,newdata=z,type='response',verbose=FALSE)) else as.numeric(predict(m$fit,newdata=z,type='response'))};write.table(as.data.frame(out,check.names=FALSE),paste0(line,'.pred'),sep='\t',row.names=FALSE,quote=FALSE);cat('BATCHDONE\n');flush.console()}
 '''.replace('LIB',str(P/'R_library')).replace('OLD',str(D)).replace('NAMES','c('+','.join(repr(n) for n in rnames)+')').replace('MODELS',str(N/'R_base'));(O/'R_service.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n);log=open(O/'R_stderr.log','w');proc=subprocess.Popen(['/usr/bin/Rscript','--vanilla',str(O/'R_service.R')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=log,text=True,env=env)
 while True:
  line=proc.stdout.readline()
  if 'READY' in line:break
  if proc.poll() is not None:raise RuntimeError('R servicefailed')
def predict(data):
 global calls
 calls+=1;pred={n:m.predict_proba(np.asarray(data))[:,1] for n,m in pymodels.items()}
 if proc:
  p=O/'batch.tsv';pd.DataFrame(data,columns=X.columns).to_csv(p,sep='\t',index=False);proc.stdin.write(str(p)+'\n');proc.stdin.flush()
  while True:
   line=proc.stdout.readline()
   if 'BATCHDONE' in line:break
   if proc.poll() is not None:raise RuntimeError('R serviceterminated')
  rd=pd.read_csv(str(p)+'.pred',sep='\t');pred.update({n:rd[n].to_numpy() for n in rnames})
 if meta is None:return pred[names[0]]
 return meta.predict_proba(pd.DataFrame({n:pred[n] for n in names}))[:,1]
ep=predict(E.to_numpy());old=pd.read_csv(N/'search/all243_external_predictions.tsv',sep='\t',index_col=0).loc[E.index,selected].to_numpy();err=float(np.max(np.abs(ep-old)));assert err<1e-8,err;np.random.seed(20261020);bg=X.sample(20,random_state=21);explainer=shap.KernelExplainer(predict,bg.to_numpy());ns=max(100,2*X.shape[1]+2);sv=np.asarray(explainer.shap_values(E.to_numpy(),nsamples=ns,l1_reg='num_features('+str(X.shape[1])+')',silent=True));pd.DataFrame(sv,index=E.index,columns=X.columns).to_csv(O/'external_gene_SHAP.tsv',sep='\t');pd.Series(np.abs(sv).mean(0),index=X.columns,name='mean_abs_SHAP').sort_values(ascending=False).to_csv(O/'gene_importance.tsv',sep='\t');add=float(np.max(np.abs(sv.sum(1)+float(explainer.expected_value)-ep)));assert add<1e-6
if proc:proc.stdin.write('STOP\n');proc.stdin.flush();proc.wait(timeout=30)
v={'status':'completed','selected':selected,'input_genes':X.shape[1],'external_samples':len(E),'nsamples':ns,'training_background':20,'prediction_max_abs_error':err,'additivity_max_abs_error':add,'method':'Actualgene-input probabilityfunction with savednativeR/Python base models andmeta; no substituteSHAPofmeta-probabilities','upstream':'Newexpandednetworkcandidatefeatures physicallyconsumed; externalpreviouslyseenexploratoryvalidation'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
