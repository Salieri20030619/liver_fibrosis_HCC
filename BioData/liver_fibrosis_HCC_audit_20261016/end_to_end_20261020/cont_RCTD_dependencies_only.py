from pathlib import Path
import subprocess,os,json
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'continuation_RCTD_dependencies_v1';L=O/'R_library'
r=r'''.libPaths(c('NEW','CLASSIC','LIB','OLD/irGSEA_dependencies_v1/R_project_library_v1','OLD/trajectory_dependencies_v1/R_project_library_v1','OLD/F6_CellChat_GSE244832_v1/R_project_library',.libPaths()))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=120,Ncpus=4)
pkgs<-c('readr','ggplot2','pals','Matrix','parallel','doParallel','foreach','quadprog','tibble','dplyr','reshape2','knitr','rmarkdown','fields','mgcv','CompQuadForm','Rfast','locfdr','metafor','data.table');missing<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)];writeLines(missing,'OUT/missing_before.txt');if(length(missing))install.packages(missing,lib='NEW');missing<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)];writeLines(missing,'OUT/missing_after.txt');if(!length(missing))writeLines('complete','OUT/DEPENDENCIES_COMPLETE.txt')
'''.replace('NEW',str(L)).replace('CLASSIC',str(P/'R_harmony_classic')).replace('LIB',str(P/'R_library')).replace('OLD',str(D)).replace('OUT',str(O));(O/'dependencies_only.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
with open(O/'dependencies_only.log','w') as f:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'dependencies_only.R')],env=env,stdout=f,stderr=subprocess.STDOUT)
print('dependencies returncode',z.returncode,flush=True)
