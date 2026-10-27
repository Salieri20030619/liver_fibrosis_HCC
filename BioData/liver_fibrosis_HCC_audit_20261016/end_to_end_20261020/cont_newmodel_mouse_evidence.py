from pathlib import Path
import json,hashlib,numpy as np,pandas as pd,anndata as ad
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'continuation_newmodel_mouse_candidates_v1';O.mkdir(exist_ok=False);mv=json.loads((P/'continuation_newnetwork_models_v2/validation.json').read_text());imp=P/('continuation_selected_model_SHAP_v1/gene_importance.tsv' if mv['status']=='reused_identical_feature_matrix' else 'continuation_newnetwork_models_v2/SHAP/gene_importance.tsv');top=pd.read_csv(imp,sep='\t',index_col=0).index.astype(str).tolist()[:3];old=P/'continuation_mouse_three_candidates_v1';prior=pd.read_csv(old/'verified_orthologs.tsv',sep='\t')
if set(top)==set(prior.human):
 v={'status':'reused_identical_top3_candidates','current_SHAP_top3':top,'importance_sha256':hashlib.sha256(imp.read_bytes()).hexdigest(),'source':str(old),'reason':'Same three verified genes and same mouse data; previous gene-expression outputs reused without recomputation. Rank order recorded separately.','confirmed_animals':None};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v));raise SystemExit(0)
M=P/'continuation_new_top3_mapping';selection=json.loads((M/'selection.json').read_text());assert selection['genes']==top;mapping=pd.read_csv(M/'new_top3_verified_orthologs.tsv',sep='\t');assert mapping.human.is_unique and set(mapping.human)<=set(top)
a=ad.read_h5ad(D/'GSE166504_S1_S7_preparation_v1/myeloid_15_30w_source_counts.h5ad');assert a.var_names.is_unique;total=np.asarray(a.X.sum(1)).ravel();z=a.obs[['source_unit','timepoint','CellType','FileName']].copy();rows=[];matched=[]
for gene in top:
 row={'human':gene}
 if gene not in set(mapping.human):row['status']='no_verified_unique_official_orthology'
 else:
  r=mapping.set_index('human').loc[gene];row.update(mouse=r.mouse,mouse_ensembl=r.mouse_ensembl)
  if r.mouse not in a.var_names:row['status']='verified_ortholog_absent_from_existing_mouse_matrix'
  else:
   values=np.asarray(a[:,r.mouse].X.toarray()).ravel();z[gene+'_mouse_log1pCP10K']=np.log1p(values/np.maximum(total,1)*1e4);matched.append(gene);row['status']='matched_and_extracted'
 rows.append(row)
pd.DataFrame(rows).to_csv(O/'candidate_mouse_mapping_status.tsv',sep='\t',index=False);mapping.to_csv(O/'verified_orthologs.tsv',sep='\t',index=False);z.to_csv(O/'all_myeloid_current_candidates.tsv.gz',sep='\t');cols=[c for c in z if c.endswith('log1pCP10K')];means=z.groupby(['source_unit','timepoint','CellType'],observed=True)[cols].mean();means.to_csv(O/'observed_source_unit_means.tsv',sep='\t');v={'status':'completed' if len(matched)==len(top) else 'partial_orthology_or_feature_coverage','current_SHAP_top3':top,'matched_candidates':matched,'importance_sha256':hashlib.sha256(imp.read_bytes()).hexdigest(),'cells':len(z),'source_units':z.source_unit.nunique(),'observed_group_rows':len(means),'confirmed_animals':None,'scope':'Existing G166504 myeloid15/30week expression with official unique mouse orthology; sourceunits not proven independent animals. No human longitudinal or incidentHCC validation claim.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
