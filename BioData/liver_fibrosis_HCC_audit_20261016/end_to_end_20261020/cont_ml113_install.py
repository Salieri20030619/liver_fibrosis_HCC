from pathlib import Path
import os,subprocess,json,datetime
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'ml113_repaired_20261020_v1';L=O/'R_library';L.mkdir(exist_ok=False);T=O/'tmp';T.mkdir();dl=O/'package_sources';dl.mkdir()
r=f""".libPaths(c('{L}','{P}/R_library','{D}/trajectory_dependencies_v1/R_project_library_v1','{D}/irGSEA_dependencies_v1/R_project_library_v1',.libPaths()))
options(timeout=180,repos=c(CRAN='https://cloud.r-project.org'))
pk<-c('e1071','randomForestSRC','gbm','xgboost')
for(p in pk){{cat('INSTALL ',p,'\\n');tryCatch(install.packages(p,lib='{L}',dependencies=NA,Ncpus=2,destdir='{dl}'),error=function(e)cat('ERROR ',conditionMessage(e),'\\n'))}}
v<-data.frame(package=pk,available=sapply(pk,requireNamespace,quietly=TRUE));v$version<-sapply(pk,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p)) else NA);write.table(v,'{O}/dependencies_after_install.tsv',sep='\\t',quote=FALSE,row.names=FALSE);writeLines(capture.output(sessionInfo()),'{O}/sessionInfo_install.txt')
"""
(O/'install_scoped.R').write_text(r);env=os.environ.copy();env['TMPDIR']=str(T);env['MAKEFLAGS']='-j2';n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
with open(O/'install.log','w') as f:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'install_scoped.R')],env=env,stdout=f,stderr=subprocess.STDOUT)
(O/'install_receipt.json').write_text(json.dumps({'R_returncode':z.returncode,'source':'https://cloud.r-project.org','destination':str(L),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
