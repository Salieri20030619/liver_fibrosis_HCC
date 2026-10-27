from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
from scipy.stats import spearmanr
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'continuation_annotation_refined_v1';q=pd.read_csv(P/'continuation_annotation_evidence_v1/complete_marker_evidence.tsv',sep='\t');q=q.drop_duplicates(['cluster','gene']).pivot(index='cluster',columns='gene',values='mean_log1pCP10k');r=ad.read_h5ad(D/'GSE136103/qc_pass_counts.h5ad');m=pd.read_csv(D/'coarse_celltype_review_v2/GSE136103_reviewed_cell_metadata.tsv.gz',sep='\t',index_col=0).loc[r.obs_names]
if 'gene_symbol_original' in r.var:r.var_names=r.var.gene_symbol_original.astype(str)
r=r[:,~r.var_names.duplicated(keep=False)].copy();genes=q.columns.intersection(r.var_names);total=np.asarray(r.X.sum(1)).ravel();x=r[:,genes].X.astype(float).multiply((1e4/np.maximum(total,1))[:,None]).tocsr();x.data=np.log1p(x.data);rows=[];means=[];names=[]
for label in sorted(m.coarse_label_reviewed.unique()):
 if label.startswith('Unresolved'):continue
 per=[]
 for donor in m.donor.unique():
  use=(m.coarse_label_reviewed.eq(label)&m.donor.eq(donor)).to_numpy()
  if use.sum()>=10:per.append(np.asarray(x[use].mean(0)).ravel())
 if per:means.append(np.mean(per,axis=0));names.append(label)
ref=pd.DataFrame(means,index=names,columns=genes);ref.to_csv(O/'external_reference_marker_profiles.tsv',sep='\t');qs=q[genes];qz=(qs-qs.mean())/qs.std().replace(0,np.nan);rz=(ref-ref.mean())/ref.std().replace(0,np.nan);good=qz.notna().all()&rz.notna().all()
for c in qs.index:
 for label in ref.index:rows.append({'cluster':int(c),'reference_label':label,'within_source_gene_centered_marker_spearman':spearmanr(qz.loc[c,good],rz.loc[label,good]).statistic,'genes':int(good.sum())})
out=pd.DataFrame(rows);out.to_csv(O/'external_reference_marker_specificity_concordance.tsv',sep='\t',index=False);(O/'reference_limitations.json').write_text(json.dumps({'status':'comparison_completed_not_ground_truth_validation','raw2000gene_profile_correlations':'Mostly negative; strong source/nuclei-vs-cell and coverage effects. Do not accept forced nearestlabel.','sensitivity':'Within-source gene-centered81canonicalmarker profiles; reference source balanced within marker-reviewed classes; retains both analyses.','missing_reference_lineages':['Hepatocyte','Stromal_fibroblast_as_distinct_label'],'reference_labels_are_own_prior_marker_reviews':True},indent=2));print(out.sort_values('within_source_gene_centered_marker_spearman',ascending=False).groupby('cluster').head(1).sort_values('cluster').to_json(orient='records'),flush=True)
