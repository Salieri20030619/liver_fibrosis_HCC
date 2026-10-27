from pathlib import Path
import os,subprocess,json,hashlib,tarfile,traceback,sys
P=Path(__file__).parent;mode=sys.argv[1];O=P/('dependency_'+mode+'_repair_v5');O.mkdir(exist_ok=False)
def run(args,env=None):
 with open(O/'execution.log','a',buffering=1) as f:return subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,env=env).returncode
try:
 if mode=='cytotrace':
  url='https://codeload.github.com/gunsagargulati/CytoTRACE/tar.gz/refs/heads/master';f=O/'CytoTRACE_master.tar.gz';ret=run(['curl','--noproxy','*','-fL','--connect-timeout','15','--max-time','1200','--max-filesize','100000000',url,'-o',str(f)]);assert ret==0,('download_failed',ret)
  with tarfile.open(f,'r:gz') as tar:members=tar.getnames();assert any(x.endswith('/DESCRIPTION') for x in members)
  (O/'source.json').write_text(json.dumps({'url':url,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}));lib=P/'trajectory_dependencies_v1/R_project_library_v1';s=f"lib<-'{lib}'\n.libPaths(c(lib,.Library))\noptions(repos=c(CRAN='https://cloud.r-project.org'),timeout=600)\nremotes::install_local('{f}',lib=lib,upgrade='never',dependencies=NA)\nstopifnot(requireNamespace('CytoTRACE',quietly=TRUE))\nwriteLines(as.character(packageVersion('CytoTRACE')),'{O/'installed_version.txt'}')\n";env=os.environ.copy()
 else:
  import time
  while Path('/proc/229805').exists():time.sleep(15)
  prefix=P/'dependency_irgsea_repair_v3/native_dependencies';env=os.environ.copy();native=prefix/'lib/libstdc++.so.6';assert native.exists();env['LD_PRELOAD']=str(native)+':'+str(prefix/'lib/libdeflate.so.0')
  env['PKG_CONFIG_PATH']=str(prefix/'lib/pkgconfig');env['PATH']=str(prefix/'bin')+os.pathsep+env['PATH'];env['LD_LIBRARY_PATH']=str(prefix/'lib')+os.pathsep+env.get('LD_LIBRARY_PATH','');env['USE_BUNDLED_LIBUV']='1';lib=P/'irGSEA_dependencies_v1/R_project_library_v1';fallback=P/'trajectory_dependencies_v1/R_project_library_v1';s=f"lib<-'{lib}'\n.libPaths(c(lib,'{fallback}',.Library))\noptions(timeout=600,Ncpus=4,repos=BiocManager::repositories())\ninstall.packages('https://cran.r-project.org/src/contrib/Archive/gghalves/gghalves_0.1.4.tar.gz',repos=NULL,type='source',lib=lib)\nBiocManager::install(c('magick','SpatialExperiment','GSVA'),lib=lib,ask=FALSE,update=FALSE,Ncpus=4)\nremotes::install_github('chuiqin/irGSEA',lib=lib,upgrade='never',dependencies=NA)\nstopifnot(requireNamespace('irGSEA',quietly=TRUE))\nwriteLines(as.character(packageVersion('irGSEA')),'{O/'installed_version.txt'}')\n"
 q=O/'install.R';q.write_text(s);ret=run(['/usr/bin/Rscript','--vanilla',str(q)],env);(O/'result.json').write_text(json.dumps({'returncode':ret,'complete':(O/'installed_version.txt').exists()}))
except Exception:(O/'error.txt').write_text(traceback.format_exc());raise
