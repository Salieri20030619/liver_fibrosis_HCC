import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
import pathlib,sys,tarfile,gzip,io,json,hashlib,datetime
import numpy as np,pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits
ROOT=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');g=sys.argv[1]
def check_matrix(blob,genes,barcodes):
 with threadpool_limits(limits=2):a=mmread(gzip.GzipFile(fileobj=io.BytesIO(blob)))
 gids=[x.split('\t')[0] for x in genes];r={'rows':a.shape[0],'columns':a.shape[1],'nnz':a.nnz,'features':len(genes),'barcodes':len(barcodes),'duplicate_feature_ids':len(gids)-len(set(gids)),'duplicate_barcodes':len(barcodes)-len(set(barcodes)),'nonfinite_values':int((~np.isfinite(a.data)).sum()),'negative_values':int((a.data<0).sum()),'noninteger_values':int((a.data!=np.floor(a.data)).sum()),'zero_library_cells':int((np.asarray(a.sum(axis=0)).ravel()==0).sum()),'matrix_gzip_sha256':hashlib.sha256(blob).hexdigest(),'gene_id_order_sha256':hashlib.sha256('\n'.join(gids).encode()).hexdigest()}
 r['structural_pass']=a.shape==(len(genes),len(barcodes)) and not any(r[k] for k in ['duplicate_feature_ids','duplicate_barcodes','nonfinite_values','negative_values','noninteger_values']);return r
rows=[]
if g=='GSE136103':
 p=ROOT/'GSE136103_RAW.tar.complete';assert p.exists();meta=pd.read_csv(ROOT/'GSE136103_library_donor_mapping.tsv',sep='\t',keep_default_na=False).set_index('GSM')
 with tarfile.open(p,'r:') as t:
  members={m.name:m for m in t.getmembers() if m.isfile()}
  for name,m in members.items():
   if not name.endswith('_matrix.mtx.gz'):continue
   pre=name[:-len('_matrix.mtx.gz')];gn=pre+'_genes.tsv.gz';gn=gn if gn in members else pre+'_features.tsv.gz';bn=pre+'_barcodes.tsv.gz'
   gb=t.extractfile(members[gn]).read();bb=t.extractfile(members[bn]).read();genes=gzip.decompress(gb).decode().splitlines();bc=gzip.decompress(bb).decode().splitlines();r=check_matrix(t.extractfile(m).read(),genes,bc);gsm=name.split('_')[0];r.update(GSM=gsm,filename=name,barcodes_sha256=hashlib.sha256(bb).hexdigest(),features_sha256=hashlib.sha256(gb).hexdigest())
   if gsm in meta.index:r.update({k:meta.loc[gsm,k] for k in ['derived_donor_id','include_human_liver','organism_ch1','source_name_ch1','cause of liver disease']})
   rows.append(r);print(json.dumps(r,default=str),flush=True)
elif g=='GSE189175':
 p=ROOT/'GSE189175_counts_matrix.mtx.gz.complete';genes=gzip.open(ROOT/'GSE189175_counts_features.tsv.gz','rt').read().splitlines();bc=gzip.open(ROOT/'GSE189175_counts_barcodes.tsv.gz','rt').read().splitlines();r=check_matrix(p.read_bytes(),genes,bc);m=pd.read_csv(ROOT/'GSE189175_meta_data.tsv.gz',sep='\t',keep_default_na=False);r.update(metadata_cells=len(m),matrix_only_barcodes=len(set(bc)-set(m.iloc[:,0])),metadata_only_barcodes=len(set(m.iloc[:,0])-set(bc)),same_order=bc==m.iloc[:,0].tolist(),patients=m.Patient.nunique(),samples=m.SampleID.nunique(),unresolved_gene_id_rows=sum(x.split('\t')[0] in ['NA','N/A',''] for x in genes));rows.append(r)
pd.DataFrame(rows).to_csv(ROOT/(g+'_full_matrix_audit.tsv'),sep='\t',index=False)
result={'resource':g,'all_matrices_structurally_pass':all(x['structural_pass'] for x in rows),'matrices':len(rows),'cells_total_in_files':sum(x['columns'] for x in rows),'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'rows':rows}
(ROOT/(g+'_full_validation.json')).write_text(json.dumps(result,indent=2,default=str));print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
