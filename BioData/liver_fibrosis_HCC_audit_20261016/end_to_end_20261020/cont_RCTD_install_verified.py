from pathlib import Path
import time,subprocess,os,json,hashlib
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'continuation_RCTD_dependencies_v1';archive=O/'spacexr_runtime_complete.tar.gz'
for i in range(180):
 if (O/'DEPENDENCIES_COMPLETE.txt').exists() and (O/'RCTD_archive_receipt.json').exists():break
 time.sleep(10)
else:raise RuntimeError('Dependencies or verifiedsource pending')
receipt=json.loads((O/'RCTD_archive_receipt.json').read_text());assert archive.stat().st_size==receipt['bytes'];assert hashlib.sha256(archive.read_bytes()).hexdigest()==receipt['sha256']
r=r'''.libPaths(c('NEW','CLASSIC','LIB','OLD/irGSEA_dependencies_v1/R_project_library_v1','OLD/trajectory_dependencies_v1/R_project_library_v1','OLD/F6_CellChat_GSE244832_v1/R_project_library',.libPaths()));install.packages('ARCHIVE',repos=NULL,type='source',lib='NEW');stopifnot(requireNamespace('spacexr',quietly=TRUE));writeLines(as.character(packageVersion('spacexr')),'OUT/COMPLETE.txt');writeLines(capture.output(sessionInfo()),'OUT/sessionInfo.txt')
'''.replace('NEW',str(O/'R_library')).replace('CLASSIC',str(P/'R_harmony_classic')).replace('LIB',str(P/'R_library')).replace('OLD',str(D)).replace('ARCHIVE',str(archive)).replace('OUT',str(O));(O/'install_verified.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
with open(O/'install_verified.log','w') as f:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'install_verified.R')],env=env,stdout=f,stderr=subprocess.STDOUT)
print('verified install returncode',z.returncode,flush=True)
