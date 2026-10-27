from pathlib import Path
import os,sys,json,time,subprocess,datetime,traceback
P=Path(__file__).parent;O=P/'F5_GSE244832_regulon_v1';V=P/'F5_setup_v1/pyscenic_venv';REF=P.parent/'reference_cisTarget_hg38_mc9nr_v3';AUX=P.parent/'reference_cisTarget_aux_v1'
env=os.environ.copy()
tmp=P/'F5_setup_v1/runtime_tmp_v2';tmp.mkdir(exist_ok=True)
for k in ['TMPDIR','TMP','TEMP']:env[k]=str(tmp)
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='1'
def state(stage,**kw):
 v={'stage':stage,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw};(O/'pipeline_state.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
def waitfor(path,errors):
 start=time.monotonic()
 while not path.exists():
  for q in errors:
   if q.exists():raise RuntimeError('Dependency failed: '+str(q)+'\n'+q.read_text()[-2500:])
  if time.monotonic()-start>14400:raise RuntimeError('Dependency still pending after4hours: '+str(path))
  time.sleep(15)
def run(name,args):
 with open(O/(name+'.log'),'w') as log:
  r=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env,cwd=O);state(name,pid=r.pid,command=args);ret=r.wait()
 if ret:raise RuntimeError(name+' failed, returncode='+str(ret))
 (O/(name+'_complete.json')).write_text(json.dumps({'returncode':ret,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':args},indent=2))
try:
 assert json.loads((O/'GRN_output_validation.json').read_text())['status']=='passed'
 assert json.loads((O/'reference_validation.json').read_text())['status']=='passed'
 refs=sorted(str(q) for q in REF.glob('*.feather'));assert len(refs)==2
 run('motif_pruning',[str(V/'bin/pyscenic'),'ctx',str(O/'adjacencies.tsv'),*refs,'--annotations_fname',str(AUX/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'),'--expression_mtx_fname',str(O/'discovery_raw_counts.tsv'),'--mode','custom_multiprocessing','--num_workers','4','--output',str(O/'regulons.csv')])
 run('AUCell',[str(V/'bin/pyscenic'),'aucell',str(O/'projection_raw_counts.tsv'),str(O/'regulons.csv'),'--num_workers','4','--seed','0','--output',str(O/'aucell_matrix.csv')])
 run('CSI_and_donor_summary',[str(V/'bin/python'),str(P/'stage1_F5_CSI_postprocess_v1.py')]);state('complete')
except Exception:(O/'pipeline_error.txt').write_text(traceback.format_exc());state('failed',error=traceback.format_exc());raise
