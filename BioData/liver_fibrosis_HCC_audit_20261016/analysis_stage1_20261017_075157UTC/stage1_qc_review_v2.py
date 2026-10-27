import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,gzip,re,datetime,numpy as np,pandas as pd,anndata as ad,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'QC_review_v2';O.mkdir(exist_ok=False)
ref=json.loads((P/'ensembl_MT_gene_reference.json').read_text());feats=pd.read_csv(P.parent/'GSE189175_counts_features.tsv.gz',sep='\t',header=None,names=['source_id','symbol'],keep_default_na=False)
feats['stable_id']=feats.source_id.str.replace(r'[.][0-9]+$','',regex=True);ids={g['id'] for g in ref};symbols={g.get('external_name','') for g in ref};hits=feats[feats.stable_id.isin(ids)|feats.symbol.isin(symbols)]
ev={'feature_rows':len(feats),'versioned_Ensembl_rows':int(feats.source_id.str.match(r'^ENSG[0-9]+[.][0-9]+$').sum()),'reference_genes':len(ids),'mitochondrial_matches':hits.to_dict('records'),'feature_columns':'versioned Ensembl gene ID, gene symbol','reference_scope':'current official human MT stable IDs and symbols; original annotation release unconfirmed','interpretation':'No missing mitochondrial counts can be reconstructed by annotation mapping when the feature rows are absent','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'GSE189175_MT_feature_audit.json').write_text(json.dumps(ev,indent=2));pd.DataFrame([{'stable_id':g['id'],'symbol':g.get('external_name'),'present':g['id'] in set(feats.stable_id)} for g in ref]).to_csv(O/'MT_reference_presence.tsv',sep='\t',index=False)
assert len(hits)==0,'MT features found: recompute metric from original counts instead of masking'
d=O/'GSE189175';d.mkdir();source=P/'GSE189175';obs=pd.read_csv(source/'cell_QC_and_source_metadata.tsv.gz',sep='\t',index_col=0,keep_default_na=False)
obs['pct_mt_recomputed']=np.nan;obs['mt_metric_available']=False;obs['mt_QC_applied']=False;obs['fail_pct_mt_recomputed']=pd.NA
obs.to_csv(d/'cell_QC_and_source_metadata.tsv.gz',sep='\t',na_rep='NA')
a=ad.read_h5ad(source/'qc_pass_counts.h5ad');a.obs['pct_mt_recomputed']=np.nan;a.obs['mt_metric_available']=False;a.obs['mt_QC_applied']=False
if 'fail_pct_mt_recomputed' in a.obs:del a.obs['fail_pct_mt_recomputed']
a.uns['mitochondrial_QC']='Unavailable: source feature table contains none of 37 official MT genes; not applied. pct_mt is NaN, never biological zero.'
a.uns['legacy_source']=str(source/'qc_pass_counts.h5ad');a.write_h5ad(d/'qc_pass_counts_MT_NA.h5ad',compression='gzip');del a
fig,axs=plt.subplots(1,3,figsize=(12,3.5))
for ax,c in zip(axs[:2],['total_counts_recomputed','n_genes_recomputed']):ax.hist(obs[c].astype(float),bins=80);ax.set_title(c);ax.set_yscale('log')
axs[2].text(.5,.5,'Mitochondrial QC unavailable\nNot applied; source features absent',ha='center',va='center');axs[2].set_axis_off();fig.tight_layout();fig.savefig(d/'F2_preparation_QC_MT_NA.png',dpi=160);plt.close(fig)
rates=[];thresholds=[]
for g in ['GSE244832','GSE136103','GSE189175','GSE212837']:
 q=P/g/'QC_by_donor_library.tsv'
 if q.exists():
  df=pd.read_csv(q,sep='\t');df.insert(0,'cohort',g);df['retention_rate']=df.pass_cells/df.input_cells;df['mt_QC_applied']=g!='GSE189175';rates.append(df)
 q=P/g/'QC_thresholds_per_library.tsv'
 if q.exists():
  df=pd.read_csv(q,sep='\t');df.insert(0,'cohort',g);df['metric_available']=True;df['filter_applied']=~df.zero_MAD_no_filter.astype(bool)
  if g=='GSE189175':
   sel=df.metric.eq('pct_mt_recomputed');df.loc[sel,['median','MAD_unscaled','lower','upper','flagged']]=np.nan;df.loc[sel,['metric_available','filter_applied']]=False
  thresholds.append(df)
pd.concat(rates).to_csv(O/'QC_retention_by_donor_library.tsv',sep='\t',index=False);pd.concat(thresholds).to_csv(O/'QC_thresholds_all_libraries.tsv',sep='\t',index=False,na_rep='NA')
rows=[['QC_multiplier','3 MAD','baseline_explicit'],['QC_grouping','per library','reconstructed; baseline grouping unspecified'],['QC_metric_scale','raw counts/features and percent mitochondrial','reconstructed; baseline scale unspecified'],['MAD_definition','median absolute deviation, unscaled','reconstructed; baseline scaling unspecified'],['QC_tails','two-sided counts/features; upper mitochondrial','reconstructed; baseline tails unspecified'],['zero_MAD','do not apply outlier threshold; flag for review','reconstructed'],['zero_library_sum','exclude zero-total cells','computational validity requirement'],['minimum_gene_fixed_cutoff','none added','no additional gate'],['GSE189175_MT','NA; filter not applied','source metric unavailable'],['HVG','3000','baseline_explicit'],['PCA_Harmony_neighbors_UMAP','workflow retained','baseline_explicit'],['implementation','Scanpy rather than Seurat4.3','reconstruction; original R package unavailable'],['PC_neighbors_resolution','30 PCs,15 neighbors,Leiden0.5','reconstructed; baseline unspecified']]
pd.DataFrame(rows,columns=['parameter','actual_choice','provenance']).to_csv(O/'baseline_and_reconstructed_parameters.tsv',sep='\t',index=False)
result={'status':'complete','GSE189175_corrected_h5ad':str(d/'qc_pass_counts_MT_NA.h5ad'),'GSE189175_retained_cells':int(obs.qc_pass.astype(str).str.lower().eq('true').sum()),'QC_selection_changed':False,'legacy_zero_metric_outputs_superseded_for_downstream':True,'cohorts_in_retention_table':[str(df.cohort.iloc[0]) for df in rates]}
(O/'result.json').write_text(json.dumps(result,indent=2));(P/'downstream_QC_input_overrides_v2.json').write_text(json.dumps({'GSE189175':result['GSE189175_corrected_h5ad'],'instruction':'Use corrected NA metric and QC figure; do not include legacy mitochondrial zeros in summary plots; retained cells unchanged'},indent=2));print(json.dumps(result),flush=True)
