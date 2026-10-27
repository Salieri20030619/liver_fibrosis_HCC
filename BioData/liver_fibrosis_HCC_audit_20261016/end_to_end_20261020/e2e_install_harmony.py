from pathlib import Path
import subprocess,json,os
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';lib=P/'R_library';lib.mkdir(exist_ok=True)
env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
r=f".libPaths(c('{lib}','{D}/irGSEA_dependencies_v1/R_project_library_v1','{D}/trajectory_dependencies_v1/R_project_library_v1',.libPaths())); options(timeout=180,download.file.method='libcurl');install.packages(c('harmony','glmnet','mboost','plsRglm'),repos='https://cloud.r-project.org',lib='{lib}',Ncpus=4);for(n in c('harmony','glmnet','mboost','plsRglm'))cat(n,requireNamespace(n,quietly=TRUE),'\\n')"
(P/'install_methods.R').write_text(r)
z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(P/'install_methods.R')],env=env,timeout=900)
print(json.dumps({'returncode':z.returncode}));raise SystemExit(z.returncode)
