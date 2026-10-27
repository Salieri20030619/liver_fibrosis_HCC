from pathlib import Path
import os,time,json,subprocess,traceback,shutil
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';I=P/'continuation_network_communication_inputs_v1/SCENIC';O=P/'continuation_SCENIC_expanded_v1';O.mkdir(exist_ok=False);V=D/'F5_setup_v1/pyscenic_venv';REF=P.parent/'reference_cisTarget_hg38_mc9nr_v3';AUX=P.parent/'reference_cisTarget_aux_v1';tmp=O/'runtime_tmp';tmp.mkdir();env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='4'
for k in ['TMPDIR','TMP','TEMP']:env[k]=str(tmp)
def run(stage,args):
 (O/'state.json').write_text(json.dumps({'stage':stage,'command':args},indent=2))
 with open(O/(stage+'.log'),'w') as f:z=subprocess.run(args,env=env,cwd=O,stdout=f,stderr=subprocess.STDOUT)
 if z.returncode:raise RuntimeError(stage+' returncode '+str(z.returncode))
 (O/(stage+'_COMPLETE.txt')).write_text('complete');print(stage,'completed',flush=True)
try:
 for i in range(240):
  if (I.parent/'validation.json').exists() and (P/'continuation_annotation_review_acceptance.json').exists():break
  time.sleep(10)
 else:raise RuntimeError('Reviewed annotation or network input not ready')
 assert shutil.disk_usage(tmp).free>30*1024**3
 run('GRNBoost2',[str(V/'bin/pyscenic'),'grn',str(I/'discovery_raw_counts.tsv'),str(D/'F5_regulatory_resource_audit_v1/allTFs_hg38.txt'),'--num_workers','4','--seed','0','--method','grnboost2','--output',str(O/'adjacencies.tsv')])
 refs=sorted(str(p) for p in REF.glob('*.feather'));assert len(refs)==2
 run('motif_pruning',[str(V/'bin/pyscenic'),'ctx',str(O/'adjacencies.tsv'),*refs,'--annotations_fname',str(AUX/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'),'--expression_mtx_fname',str(I/'discovery_raw_counts.tsv'),'--mode','custom_multiprocessing','--num_workers','4','--output',str(O/'regulons.csv')])
 run('AUCell',[str(V/'bin/pyscenic'),'aucell',str(I/'projection_raw_counts.tsv'),str(O/'regulons.csv'),'--num_workers','4','--seed','0','--output',str(O/'aucell_matrix.csv')])
 (O/'state.json').write_text(json.dumps({'stage':'GRN_motifs_AUCell_completed_postprocess_pending'}))
except Exception:(O/'error.txt').write_text(traceback.format_exc());raise
