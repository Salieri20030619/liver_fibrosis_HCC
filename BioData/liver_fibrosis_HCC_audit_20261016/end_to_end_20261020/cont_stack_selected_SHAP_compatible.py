from pathlib import Path
import os,sys,json,pickle,subprocess
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));sys.path.insert(0,str(P/'shap_compatible'));D=P.parent/'analysis_stage1_20261017_075157UTC';I=P/'F8_F9_real_full_cohort_v1';O=P/'continuation_selected_model_SHAP_v1';O.mkdir(exist_ok=True)
import numpy as np,pandas as pd,shap
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
X=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0);E=pd.read_csv(I/'external_features.tsv',sep='\t',index_col=0);y=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0).loc[X.index,'y'];selected=json.loads((P/'continuation_stack_search_v1/validation.json').read_text())['selected']['model'];assert selected=='Stack_Lasso_0.1__glmBoost_200'
with open(P/'continuation_stack_search_v1/meta_models.pkl','rb') as f:meta=pickle.load(f)[selected]
lasso=Pipeline([('scale',StandardScaler()),('RFE',RFE(LogisticRegression(C=1,solver='liblinear',max_iter=2000),n_features_to_select=14,step=1)),('model',LogisticRegression(penalty='l1',C=.1,solver='saga',max_iter=5000,random_state=3))]).fit(X.to_numpy(),y)
r=r'''.libPaths(c('LIB','OLD/irGSEA_dependencies_v1/R_project_library_v1','OLD/trajectory_dependencies_v1/R_project_library_v1',.libPaths()));library(mboost)
model<-readRDS('MODEL');con<-file('stdin',open='r');cat('READY\n');flush.console()
repeat{line<-readLines(con,n=1);if(!length(line)||line=='STOP')break;x<-as.matrix(read.delim(line,check.names=FALSE));z<-scale(x,model$preprocessing$mean,model$preprocessing$sd);z<-as.data.frame(z[,model$preprocessing$features,drop=FALSE]);p<-as.numeric(predict(model$fit,newdata=z,type='response'));write.table(p,paste0(line,'.pred'),row.names=FALSE,col.names=FALSE,quote=FALSE);cat('BATCHDONE\n');flush.console()}
'''.replace('LIB',str(P/'R_library')).replace('OLD',str(D)).replace('MODEL',str(P/'F8_remaining_R_families/glmBoost_200.rds'))
(O/'prediction_service.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n);log=open(O/'R_service_stderr.log','w');proc=subprocess.Popen(['/usr/bin/Rscript','--vanilla',str(O/'prediction_service.R')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=log,text=True,env=env)
while True:
 line=proc.stdout.readline()
 if 'READY' in line:break
 if proc.poll() is not None:raise RuntimeError('R prediction service failed')
calls=0
def predict(data):
 global calls
 calls+=1;p=O/'prediction_batch.tsv';pd.DataFrame(data,columns=X.columns).to_csv(p,sep='\t',index=False);proc.stdin.write(str(p)+'\n');proc.stdin.flush()
 while True:
  line=proc.stdout.readline()
  if 'BATCHDONE' in line:break
  if proc.poll() is not None:raise RuntimeError('R prediction service terminated')
 rp=np.loadtxt(str(p)+'.pred');lp=lasso.predict_proba(np.asarray(data))[:,1];return meta.predict_proba(pd.DataFrame({'Lasso_0.1':lp,'glmBoost_200':np.atleast_1d(rp)}))[:,1]
ep=predict(E.to_numpy());old=pd.read_csv(P/'continuation_stack_search_v1/all243_external_predictions.tsv',sep='\t',index_col=0).loc[E.index,selected].to_numpy();err=float(np.max(np.abs(ep-old)));assert err<1e-8,err
np.random.seed(20261020);bg=X.sample(20,random_state=21);explainer=shap.KernelExplainer(predict,bg.to_numpy());values=explainer.shap_values(E.to_numpy(),nsamples=100,l1_reg='num_features(28)');values=np.asarray(values);pd.DataFrame(values,index=E.index,columns=X.columns).to_csv(O/'external112_gene_SHAP.tsv',sep='\t');pd.Series(np.abs(values).mean(0),index=X.columns,name='mean_abs_SHAP').sort_values(ascending=False).to_csv(O/'gene_importance.tsv',sep='\t');add=float(np.max(np.abs(values.sum(1)+float(explainer.expected_value)-ep)))
with open(O/'python_lasso_and_meta.pkl','wb') as f:pickle.dump({'lasso':lasso,'meta':meta},f)
proc.stdin.write('STOP\n');proc.stdin.flush();proc.wait(timeout=30);v={'status':'completed','model':selected,'cells_or_samples':'112 external bulk samples','features':28,'method':'KernelSHAP on actual28gene-to-probability function with native saved R glmBoost200 and reconstructed identical Python Lasso. Not SHAP of only2meta-probability inputs.','nsamples':100,'background_training_samples':20,'existing_prediction_max_abs_error':err,'additivity_max_abs_error':add,'prediction_batches':calls,'explainer_expected_value':float(explainer.expected_value)};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
