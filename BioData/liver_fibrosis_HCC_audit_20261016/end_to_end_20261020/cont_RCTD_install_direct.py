from pathlib import Path
import subprocess,os,json,urllib.request,hashlib,datetime
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'continuation_RCTD_dependencies_v1';O.mkdir(exist_ok=True);L=O/'R_library';L.mkdir(exist_ok=True)
url='https://codeload.github.com/dmcable/spacexr/tar.gz/refs/heads/master';archive=O/'spacexr_master.tar.gz'
try:
 z=subprocess.run(['curl','--noproxy','*','-L','--fail','--connect-timeout','15','--max-time','120','--max-filesize','100000000','-o',str(archive),url],capture_output=True,text=True,timeout=130)
 (O/'direct_download.log').write_text(z.stdout+z.stderr)
 assert z.returncode==0,z.stderr
 (O/'source_receipt.json').write_text(json.dumps({'url':url,'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
except Exception as e:(O/'source_error.txt').write_text(repr(e));raise
r=r'''.libPaths(c('NEW','CLASSIC','LIB','OLD/irGSEA_dependencies_v1/R_project_library_v1','OLD/trajectory_dependencies_v1/R_project_library_v1','OLD/F6_CellChat_GSE244832_v1/R_project_library',.libPaths()))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=180,Ncpus=4)
pkgs<-c('readr','ggplot2','pals','Matrix','parallel','doParallel','foreach','quadprog','tibble','dplyr','reshape2','knitr','rmarkdown','fields','mgcv','CompQuadForm','Rfast','locfdr','metafor','data.table')
missing<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)];writeLines(missing,'OUT/missing_before.txt');if(length(missing))install.packages(missing,lib='NEW')
install.packages('ARCHIVE',repos=NULL,type='source',lib='NEW');stopifnot(requireNamespace('spacexr',quietly=TRUE));writeLines(as.character(packageVersion('spacexr')),'OUT/COMPLETE.txt');writeLines(capture.output(sessionInfo()),'OUT/sessionInfo.txt')
'''.replace('NEW',str(L)).replace('CLASSIC',str(P/'R_harmony_classic')).replace('LIB',str(P/'R_library')).replace('OLD',str(D)).replace('OUT',str(O)).replace('ARCHIVE',str(archive))
(O/'install.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
with open(O/'install.log','w') as f:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'install.R')],env=env,stdout=f,stderr=subprocess.STDOUT)
print('install returncode',z.returncode,flush=True)
