from pathlib import Path
import os,subprocess,json,hashlib
import pandas as pd
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'ml113_repaired_20261020_v1';O.mkdir(exist_ok=False)
F=P/'continuation_newnetwork_models_v2/features'
v={'input_files':[]}
for q in F.iterdir():
 if q.is_file():v['input_files'].append({'name':q.name,'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest()})
for n in ['train_features.tsv','external_features.tsv','train_design.tsv','external_design.tsv','outer_cv_folds.tsv']:
 x=pd.read_csv(F/n,sep='\t',index_col=0);v[n]={'shape':list(x.shape),'columns':x.columns.tolist(),'head':x.head(2).to_dict()}
r=f""".libPaths(c('{P}/R_library','{D}/trajectory_dependencies_v1/R_project_library_v1','{D}/irGSEA_dependencies_v1/R_project_library_v1',.libPaths()))
pkgs<-c('glmnet','e1071','MASS','mboost','plsRglm','randomForestSRC','gbm','xgboost','pROC','jsonlite','caret')
x<-data.frame(package=pkgs,available=sapply(pkgs,requireNamespace,quietly=TRUE));x$version<-sapply(pkgs,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p)) else NA);write.table(x,'{O}/dependencies.tsv',sep='\\t',quote=FALSE,row.names=FALSE);writeLines(capture.output(sessionInfo()),'{O}/sessionInfo_inventory.txt')
"""
(O/'inventory.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'inventory.R')],env=env,capture_output=True,text=True);v.update(returncode=z.returncode,stdout=z.stdout,stderr=z.stderr);(O/'input_inventory.json').write_text(json.dumps(v,indent=2,default=str));print(json.dumps(v,default=str));print((O/'dependencies.tsv').read_text() if (O/'dependencies.tsv').exists() else '')
