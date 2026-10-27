from pathlib import Path
import subprocess,os,json
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';lib=P/'R_monocle_legacy';lib.mkdir(exist_ok=False)
env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
r=f".libPaths(c('{lib}','{P}/R_library','{D}/trajectory_dependencies_v1/R_project_library_v1','{D}/irGSEA_dependencies_v1/R_project_library_v1',.libPaths())); options(timeout=180);install.packages('https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.1.4.tar.gz',repos=NULL,type='source',lib='{lib}');stopifnot(as.character(packageVersion('dplyr'))=='1.1.4');writeLines(capture.output(sessionInfo()),'{lib}/sessionInfo.txt')"
(P/'install_monocle_legacy_dplyr.R').write_text(r)
with open(P/'monocle_legacy_dplyr_install.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla','-e',r],env=env,stdout=log,stderr=subprocess.STDOUT,timeout=500)
print(json.dumps({'returncode':z.returncode}));assert z.returncode==0
