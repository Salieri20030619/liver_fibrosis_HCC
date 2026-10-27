import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
import pathlib,tarfile,gzip,io,json,hashlib,datetime
import numpy as np,pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');p=r/'GSE212837_RAW.tar.complete';assert p.exists() and (r/'GSE212837_resume_v2/complete.json').exists()
meta=pd.read_csv(r/'GSE212837_library_donor_mapping.tsv',sep='\t',keep_default_na=False).set_index('GSM');rows=[];files=[]
with tarfile.open(p,'r:') as t:
 members={m.name:m for m in t.getmembers() if m.isfile()}
 for name,m in members.items():
  if not name.endswith('_matrix.mtx.gz'):continue
  pre=name[:-len('_matrix.mtx.gz')];gn=pre+'_features.tsv.gz';bn=pre+'_barcodes.tsv.gz';blobs={n:t.extractfile(members[n]).read() for n in [name,gn,bn]}
  for n,b in blobs.items():files.append({'file':n,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
  genes=gzip.decompress(blobs[gn]).decode().splitlines();bc=gzip.decompress(blobs[bn]).decode().splitlines();ids=[x.split('\t')[0] for x in genes]
  with threadpool_limits(limits=2):a=mmread(gzip.GzipFile(fileobj=io.BytesIO(blobs[name])))
  gsm=name.split('_')[0];v={'GSM':gsm,'file':name,'rows':a.shape[0],'columns':a.shape[1],'nnz':a.nnz,'features':len(genes),'barcodes':len(bc),'duplicate_feature_ids':len(ids)-len(set(ids)),'duplicate_barcodes':len(bc)-len(set(bc)),'unresolved_gene_ids':sum(x in ['', 'NA','N/A'] for x in ids),'nonfinite':int((~np.isfinite(a.data)).sum()),'negative':int((a.data<0).sum()),'noninteger':int((a.data!=np.floor(a.data)).sum()),'zero_total_cells':int((np.asarray(a.sum(axis=0)).ravel()==0).sum()),'feature_id_order_sha256':hashlib.sha256('\n'.join(ids).encode()).hexdigest(),'matrix_gzip_sha256':hashlib.sha256(blobs[name]).hexdigest(),'metadata_match':gsm in meta.index}
  if gsm in meta.index:v.update(meta.loc[gsm].to_dict())
  v['structural_pass']=a.shape==(len(genes),len(bc)) and v['metadata_match'] and not any(v[k] for k in ['duplicate_feature_ids','duplicate_barcodes','nonfinite','negative','noninteger'])
  rows.append(v);print(json.dumps(v),flush=True);del a,blobs
pd.DataFrame(rows).to_csv(r/'GSE212837_full_matrix_audit.tsv',sep='\t',index=False);pd.DataFrame(files).to_csv(r/'GSE212837_archive_member_manifest.tsv',sep='\t',index=False)
human=[x for x in rows if str(x['include_human_liver']).lower()=='true'];v={'resource':'GSE212837','finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_matrices_structurally_pass':all(x['structural_pass'] for x in rows),'libraries':len(rows),'human_libraries':len(human),'human_cells':sum(x['columns'] for x in human),'human_donors':len(set(x['derived_donor_id'] for x in human)),'excluded_mouse_libraries':len(rows)-len(human),'matrix_GSM_without_metadata':sorted(set(x['GSM'] for x in rows)-set(meta.index)),'metadata_GSM_without_matrix':sorted(set(meta.index)-set(x['GSM'] for x in rows)),'human_feature_order_hashes':sorted(set(x['feature_id_order_sha256'] for x in human)),'celltype_annotation_available_in_RAW_tar':False,'limitations':['RAW.tar counts/barcodes/features does not contain cell type labels from separate h5ad','Cell IDs must be prefixed by GSM before combining libraries','Repeated libraries are not independent donors','GSM6556454 title R1C2 conflicts with description/file R1C1; donor unchanged, preserve provenance'],'rows':rows}
(r/'GSE212837_full_validation.json').write_text(json.dumps(v,indent=2));print(json.dumps({k:x for k,x in v.items() if k!='rows'}),flush=True)
