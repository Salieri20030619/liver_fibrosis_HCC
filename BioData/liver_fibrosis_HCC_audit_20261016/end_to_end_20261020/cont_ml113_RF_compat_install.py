from pathlib import Path
import os,subprocess,json,hashlib,datetime,requests
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'ml113_repaired_20261020_v1';L=O/'R_randomForestSRC_3.3.2';L.mkdir(exist_ok=False)
url='https://cran.r-project.org/src/contrib/Archive/randomForestSRC/randomForestSRC_3.3.2.tar.gz';f=O/'package_sources/randomForestSRC_3.3.2.tar.gz';v={'url':url,'destination':str(L),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
try:
 s=requests.Session();s.trust_env=False;r=s.get(url,timeout=120);v['HTTP_status']=r.status_code;r.raise_for_status();f.write_bytes(r.content);v.update(bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())
 code=f""".libPaths(c('{L}','{O}/R_library','{P}/R_library','{D}/trajectory_dependencies_v1/R_project_library_v1','{D}/irGSEA_dependencies_v1/R_project_library_v1',.libPaths()));install.packages('{f}',lib='{L}',repos=NULL,type='source');stopifnot(requireNamespace('randomForestSRC',quietly=TRUE),'var.select'%in%getNamespaceExports('randomForestSRC'));writeLines(as.character(packageVersion('randomForestSRC')),'{O}/RF_compat_version.txt')"""
 (O/'install_RF_compat.R').write_text(code);env=os.environ.copy();env['TMPDIR']=str(O/'tmp');env['MAKEFLAGS']='-j2';n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
 with open(O/'install_RF_compat.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'install_RF_compat.R')],env=env,stdout=log,stderr=subprocess.STDOUT)
 v['returncode']=z.returncode;v['exported_var_select_verified']=(O/'RF_compat_version.txt').exists()
except Exception as e:v['error']=repr(e)
(O/'RF_compat_install_receipt.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
