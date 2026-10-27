from pathlib import Path
import os,time,gzip,json,datetime,numpy as np,pandas as pd,traceback
P=Path(__file__).parent;O=P/'gap_download';f=O/'GSE166504_cell_raw_counts.20220204.txt.gz'
try:
 deadline=time.time()+21600
 while not f.exists():
  if time.time()>deadline:raise TimeoutError('Counts download still incomplete after 6h; validation not performed')
  try:os.kill(3432264,0)
  except ProcessLookupError:raise RuntimeError('Download process exited without complete counts; inspect gap_download.log and audit.json')
  time.sleep(45)
 m=pd.read_csv(O/'GSE166504_cell_metadata.20220204.tsv.gz',sep='\t');matched=[];genes=[];sums=np.zeros(len(m));n=0
 with gzip.open(f,'rt') as inp:
  header=inp.readline().rstrip('\r\n').split('\t');cols=header[1:] if len(header)==len(m)+1 else header
  assert len(cols)==len(m),(len(cols),len(m))
  for sep in ['','_',':','|','-','.']:
   for order in [('FileName','CellID'),('CellID','FileName')]:
    expected=(m[order[0]].astype(str)+sep+m[order[1]].astype(str)).tolist()
    if cols==expected:matched.append({'order':order,'separator':sep})
  assert matched,{'column_alignment_failed':True,'matrix_first3':cols[:3],'metadata_first3':m.head(3).to_dict('records')}
  for line in inp:
   gene,values=line.rstrip('\r\n').split('\t',1);a=np.fromstring(values,sep='\t');assert a.size==len(m),(n,a.size);assert np.isfinite(a).all() and (a>=0).all(),n;sums+=a;genes.append(gene);n+=1
 out={'genes':n,'cells':len(cols),'unique_gene_ids':len(set(genes)),'unique_cell_ids':len(set(cols)),'metadata_exact_column_alignment':matched,'nonnegative_finite_counts':True,'positive_library_cells':int((sums>0).sum()),'target_gene_presence':{g:g in genes for g in ['Runx2','Plg','Pard3']},'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'unit_check':'see GSE166504_metadata_unit_check.json; libraries not independent animals'}
 (O/'GSE166504_counts_metadata_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
except Exception:(O/'GSE166504_validation_error.txt').write_text(traceback.format_exc());raise
