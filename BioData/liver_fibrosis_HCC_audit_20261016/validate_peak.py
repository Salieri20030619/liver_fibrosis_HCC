import pathlib,gzip,json,datetime,itertools,hashlib,numpy as np
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');p=r/'GSE227265_BinarizedPeaks.mtx.gz.complete'
assert p.exists(),'Full verified file absent'
with gzip.open(p,'rt') as f:
 header=f.readline().strip();line=f.readline()
 while line.startswith('%'):line=f.readline()
 nr,nc,nnz=map(int,line.split());count=0;invalid=0;order_breaks=0;duplicates=0;last=-1;counts=np.zeros(nc,dtype=np.int64)
 for batch in iter(lambda:list(itertools.islice(f,500000)),[]):
  a=np.fromstring(''.join(batch),sep=' ',dtype=np.int64).reshape(-1,2);assert len(a)==len(batch)
  invalid+=int(((a[:,0]<1)|(a[:,0]>nr)|(a[:,1]<1)|(a[:,1]>nc)).sum());keys=(a[:,1]-1)*nr+a[:,0]-1
  order_breaks+=int((np.diff(keys)<0).sum())+int(keys[0]<last);duplicates+=int((np.diff(keys)==0).sum())+int(keys[0]==last);last=int(keys[-1]);count+=len(a)
  counts+=np.bincount(a[:,1]-1,minlength=nc)
v={'resource':'GSE227265','file':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'header':header,'rows':nr,'columns':nc,'declared_entries':nnz,'parsed_entries':count,'invalid_indices':invalid,'column_major_order_breaks':order_breaks,'adjacent_duplicate_coordinates':duplicates,'duplicate_check_complete':order_breaks==0,'empty_columns':int((counts==0).sum()),'structural_pass':count==nnz and invalid==0 and order_breaks==0 and duplicates==0,'semantic_alignment_pass':False,'missing_fields':['peak genomic coordinates per row','cell IDs in matrix column order','cell to donor / HCC versus CCA map'],'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(r/'GSE227265_peak_full_validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
