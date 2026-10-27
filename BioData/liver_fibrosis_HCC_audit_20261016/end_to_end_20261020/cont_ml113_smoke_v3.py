from pathlib import Path
import os,subprocess,json,time,datetime,hashlib
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'ml113_repaired_20261020_v1'
assert (O/'smoke_receipt.json').exists() and not (O/'full_launch.json').exists()
for _ in range(360):
 if (O/'RF_compat_install_receipt.json').exists():break
 time.sleep(5)
else:raise RuntimeError('RF compatibility installation not terminal')
r=json.loads((O/'RF_compat_install_receipt.json').read_text());assert r.get('exported_var_select_verified'),r
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
with open(O/'smoke_v3.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'code/worker.R'),str(O),'smoke_v3'],env=env,stdout=log,stderr=subprocess.STDOUT)
v={'returncode':z.returncode,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (O/'code').glob('*.R')}}
(O/'smoke_v3_receipt.json').write_text(json.dumps(v,indent=2))
if z.returncode==0:
 a=json.loads((O/'jobs/smoke_v3/validation.json').read_text())
 if a['completed']==a['attempted']:(O/'smoke_acceptance.json').write_text(json.dumps({'job':'smoke_v3','receipt':v,'validation':a},indent=2))
print(json.dumps(v),flush=True)
