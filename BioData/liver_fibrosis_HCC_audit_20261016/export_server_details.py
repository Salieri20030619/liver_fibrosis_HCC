import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
import pathlib,gzip,csv,json,datetime
import pandas as pd,numpy as np
from scipy.io import mmread
from threadpoolctl import threadpool_limits
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');out=r/'review_details';out.mkdir(exist_ok=False)
allcols=[]
for slide,gsm in [(1,'GSM8855042'),(2,'GSM8855043')]:
 with gzip.open(r/f'{gsm}_Slide{slide}_exprMat_file.csv.gz','rt') as f:cols=next(csv.reader(f))[2:]
 allcols.append(cols)
assert allcols[0]==allcols[1]
rows=[]
for i,n in enumerate(allcols[0],3):
 typ='system_control' if n.startswith('SystemControl') else ('negative_probe' if 'NegPrb' in n or 'Negative' in n else 'gene_probe')
 rows.append({'csv_column_1based':i,'probe_name':n,'category':typ,'multigene_label':'/' in n,'one_to_one_gene_assumed':False if '/' in n or typ!='gene_probe' else 'not independently validated against vendor panel','source_slides':'GSM8855042;GSM8855043'})
pd.DataFrame(rows).to_csv(out/'GSE292268_1207_probe_dictionary.tsv',sep='\t',index=False)
features=gzip.open(r/'GSE189175_counts_features.tsv.gz','rt').read().splitlines();ix=[i for i,s in enumerate(features) if s.split('\t')[0] in ['NA','','N/A']];assert ix==[57536]
bc=gzip.open(r/'GSE189175_counts_barcodes.tsv.gz','rt').read().splitlines();m=pd.read_csv(r/'GSE189175_meta_data.tsv.gz',sep='\t',keep_default_na=False).set_index('Cell') if False else pd.read_csv(r/'GSE189175_meta_data.tsv.gz',sep='\t',keep_default_na=False)
assert bc==m.iloc[:,0].tolist()
with threadpool_limits(limits=2),gzip.open(r/'GSE189175_counts_matrix.mtx.gz.complete','rb') as f:a=mmread(f).tocsr()
v=a[ix[0],:].toarray().ravel();stats=[]
for sid,b in m.groupby('SampleID'):
 x=v[b.index.to_numpy()];stats.append({'SampleID':sid,'Patient':b.Patient.iloc[0],'Tissue':b.Tissue.iloc[0],'cells':len(x),'NA_feature_count_sum':float(x.sum()),'cells_with_NA_feature':int((x>0).sum()),'NA_feature_count_max':float(x.max())})
pd.DataFrame(stats).to_csv(out/'GSE189175_NA_feature_counts_by_tissue.tsv',sep='\t',index=False)
(out/'GSE189175_NA_feature_detail.json').write_text(json.dumps({'feature_row_1based':ix[0]+1,'original_feature_line':features[ix[0]],'matrix_row_index_0based':ix[0],'matrix_shape':a.shape,'row_total_counts':float(v.sum()),'nonzero_cells':int((v>0).sum()),'max_count':float(v.max()),'all_zero':bool((v==0).all()),'decision':'retain flagged unresolved feature; no unverified gene-symbol substitution','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
print(json.dumps({'probe_dictionary_rows':len(rows),'multi_gene_labels':sum(x['multigene_label'] for x in rows),'NA_feature_total':float(v.sum()),'NA_feature_nonzero_cells':int((v>0).sum())}),flush=True)
