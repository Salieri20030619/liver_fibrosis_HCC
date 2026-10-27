from pathlib import Path
import time,json,hashlib,numpy as np,pandas as pd,scanpy as sc
from scipy.stats import spearmanr
P=Path(__file__).parent;O=P/'continuation_single_slice_candidate_projection_v1';O.mkdir(exist_ok=False);C=P/'continuation_SCENIC_candidates_v1';L=P/'continuation_ligand_to_GSEA_v1';R=P/'continuation_RCTD_verified_export_v3'
for i in range(1440):
 if all(p.exists() for p in [C/'validation.json',L/'validation.json',R/'validation.json']):break
 if (P/'continuation_SCENIC_expanded_v1/error.txt').exists():raise RuntimeError('SCENIC upstream failed')
 time.sleep(10)
else:raise RuntimeError('Upstream pending4hours')
ligand=json.loads((L/'validation.json').read_text()).get('ligand');genes=sorted(set(pd.read_csv(C/'F5G_candidates.tsv',sep='\t').TF)|({ligand} if ligand else set()))
f=next((P/'continuation_RCTD_input').glob('*.h5'));a=sc.read_10x_h5(f);w=pd.read_csv(R/'weights_row_normalized.tsv',sep='\t',index_col=0);co=pd.read_csv(R/'spot_coordinates_and_proportions.tsv',sep='\t',index_col=0);assert w.index.isin(a.obs_names).all();a=a[w.index].copy();total=np.asarray(a.X.sum(axis=1)).ravel();expr={};mapping=[]
for g in genes:
 ids=np.flatnonzero(a.var_names==g);mapping.append({'gene':g,'matching_features':len(ids),'status':'unique' if len(ids)==1 else 'absent' if len(ids)==0 else 'ambiguous_excluded'})
 if len(ids)==1:expr[g]=np.log1p(np.asarray(a.X[:,ids].toarray()).ravel()/np.maximum(total,1)*1e4)
e=pd.DataFrame(expr,index=w.index);e.to_csv(O/'candidate_spot_log1pCP10K.tsv',sep='\t');co.loc[w.index].to_csv(O/'aligned_spot_coordinates_RCTD_proportions.tsv',sep='\t');pd.DataFrame(mapping).to_csv(O/'gene_mapping.tsv',sep='\t',index=False);rows=[]
for gene in e:
 for celltype in w:
  rho,_=spearmanr(e[gene],w[celltype]);rows.append({'feature1':gene,'feature2':celltype,'kind':'expression_RCTD_fraction','spots':len(w),'rho_descriptive':rho})
for i,t in enumerate(w):
 for u in w.columns[i+1:]:
  rho,_=spearmanr(w[t],w[u]);rows.append({'feature1':t,'feature2':u,'kind':'RCTD_fraction_fraction','spots':len(w),'rho_descriptive':rho})
if ligand in e:
 for gene in e:
  if gene!=ligand:
   rho,_=spearmanr(e[gene],e[ligand]);rows.append({'feature1':gene,'feature2':ligand,'kind':'TF_ligand_expression','spots':len(w),'rho_descriptive':rho})
pd.DataFrame(rows).to_csv(O/'single_slice_spot_cooccurrence_descriptive.tsv',sep='\t',index=False)
v={'status':'completed_single_slice_descriptive_projection','spots':len(w),'genes_requested':len(genes),'genes_unique_matched':len(e.columns),'ligand':ligand,'source_candidate_sha256':hashlib.sha256((C/'F5G_candidates.tsv').read_bytes()).hexdigest(),'scope':'One existing GSE203612 LIHC section; final55record7broadtype reference; numeric candidate expression and RCTD coordinates only, no formal plots. Correlations have no iidspot significance test; adjacent spots are dependent; fractions compositional.','limits':'Not MASH etiology evidence, not seven section replication, not cell-specific ligand expression, no causal co-localization or human longitudinal progression claim. Fibrosis-to-HCC reference domain shift remains.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
