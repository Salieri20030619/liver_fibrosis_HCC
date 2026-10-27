from pathlib import Path
import subprocess,os,json,time,datetime,hashlib
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'ml113_repaired_20261020_v1'
acceptance=json.loads((O/'smoke_acceptance.json').read_text());sm=acceptance['validation'];assert sm['metamorphic_preprocessing'] and sm['ID_order_rejection'] and sm['fixed_AUC_direction'] and sm['forward_scope_test']
assert all(hashlib.sha256((O/'code'/n).read_bytes()).hexdigest()==h for n,h in acceptance['receipt']['code_sha256'].items())
assert sm['attempted']==sm['completed'], 'Resolve smoke implementation/dependency errors before fullrun'
assert not (O/'full_launch.json').exists()
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1');n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
jobs=['outer_'+str(k) for k in range(10)]+['full'];active={};finished=[];launched=[]
receipt={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'max_parallel_R_jobs':3,'code_sha256':{q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in (O/'code').glob('*.R')},'jobs':launched,'finished':finished}
(O/'full_launch.json').write_text(json.dumps(receipt,indent=2))
while jobs or active:
 while jobs and len(active)<3:
  name=jobs.pop(0);assert not (O/'jobs'/name).exists()
  log=open(O/(name+'.log'),'w')
  z=subprocess.Popen(['/usr/bin/Rscript','--vanilla',str(O/'code/worker.R'),str(O),name],env=env,stdout=log,stderr=subprocess.STDOUT)
  active[name]=(z,log);launched.append({'job':name,'pid':z.pid,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()});print('START',name,z.pid,flush=True)
 for name,(z,log) in list(active.items()):
  if z.poll() is not None:
   log.close();finished.append({'job':name,'returncode':z.returncode,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()});del active[name];print('END',name,z.returncode,flush=True)
 (O/'full_launch.json').write_text(json.dumps(receipt,indent=2))
 if active:time.sleep(10)
v={'status':'all_workers_terminated','jobs':finished,'all_returncodes_zero':all(x['returncode']==0 for x in finished),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'workers_terminal.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
