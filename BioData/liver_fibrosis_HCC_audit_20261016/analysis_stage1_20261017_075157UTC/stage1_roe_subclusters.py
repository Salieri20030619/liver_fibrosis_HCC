import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,datetime,traceback,numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'F2_Roe_F3_F4_subclusters_v1';O.mkdir(exist_ok=False)
def log(task,status,**kw):
 x={'task':task,'status':status,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw};(O/'state.json').write_text(json.dumps(x));print(json.dumps(x),flush=True)
manifest={'baseline_panels':['F2E/F: observed/expected ratios','F3A/B and F4A/B: independent mesenchymal/myeloid subclustering and markers'],'preserved_methods':['Ro/e=observed/(row total*column total/grand total)','top3000 HVG, PCA, Harmony, graph clustering, UMAP'],'reconstruction':['Scanpy instead of Seurat4.3','30PC,15neighbors,Leiden0.5,UMAPmin_dist0.5,seed0'],'disease_scope':'Within study; GSE244832 F0/F3/F4 after D9/D10 exclusion; GSE136103 healthy/fibrotic separately within each CD45 sorting stratum','statistics':'Ro/e descriptive, not donor-level significance; original Wilcoxon rank-sum retained for future applicable two-group tests; no new significance gate','subcluster_scope':'Mesenchymal includes stromal/mural without assuming HSC subtype; no DTNA or RUNX2 axis assumption; no forced number of subclusters'}
(O/'analysis_manifest.json').write_text(json.dumps(manifest,indent=2))
for g in ['GSE244832','GSE136103']:
 obs=pd.read_csv(P/'coarse_celltype_review_v2'/(g+'_reviewed_cell_metadata.tsv.gz'),sep='\t',index_col=0,keep_default_na=False)
 eligible=obs.disease_comparison_eligible.astype(str).str.lower().eq('true');clinical=obs.loc[eligible].copy();clinical['comparison_group']=('F'+clinical.F_stage_from_Figure4C.astype(str)) if g=='GSE244832' else clinical['disease status']
 for sort,sub in clinical.groupby('sorting',observed=True):
  tab=pd.crosstab(sub.coarse_label_reviewed,sub.comparison_group);expected=np.outer(tab.sum(axis=1),tab.sum(axis=0))/tab.to_numpy().sum();roe=tab/expected;key=g+'_'+sort.replace('+','positive').replace('-','negative');tab.to_csv(O/(key+'_observed.tsv'),sep='\t');roe.to_csv(O/(key+'_RoE.tsv'),sep='\t');sub.groupby(['donor','comparison_group','coarse_label_reviewed']).size().rename('cells').reset_index().to_csv(O/(key+'_donor_cell_counts.tsv'),sep='\t',index=False)
  fig,ax=plt.subplots(figsize=(7,max(4,len(tab)*.3)));im=ax.imshow(roe.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=0,vmax=2);ax.set_xticks(range(len(tab.columns)),tab.columns);ax.set_yticks(range(len(tab.index)),tab.index);ax.set_title(key+' | descriptive Ro/e; not significance');fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(O/(key+'_F2_RoE.png'),dpi=160);plt.close(fig)
 for compartment in ['Mesenchymal','Myeloid']:
  task=g+'_'+compartment
  try:
   d=O/task;d.mkdir();mask=obs.coarse_label_reviewed.str.startswith(compartment);src=ad.read_h5ad(P/g/'qc_pass_counts.h5ad',backed='r');assert obs.index.tolist()==src.obs_names.tolist();a=src[mask.to_numpy(),:].to_memory();src.file.close();a.obs=obs.loc[mask].copy();log(task,'loaded',cells=a.n_obs,donors=a.obs.donor.nunique())
   sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);a.raw=a.copy();sc.pp.highly_variable_genes(a,n_top_genes=3000,flavor='seurat');a.var.to_csv(d/'HVG_selection.tsv',sep='\t');a=a[:,a.var.highly_variable].copy();sc.pp.scale(a,max_value=10);sc.tl.pca(a,n_comps=30,svd_solver='arpack',random_state=0);np.save(d/'PCA_unadjusted.npy',a.obsm['X_pca'])
   z=np.asarray(hm.run_harmony(a.obsm['X_pca'],a.obs,'library',random_state=0,max_iter_harmony=20).Z_corr)
   if z.shape!=a.obsm['X_pca'].shape:z=z.T
   assert z.shape==a.obsm['X_pca'].shape;a.obsm['X_pca_harmony']=z;np.save(d/'PCA_harmony.npy',z);sc.pp.neighbors(a,n_neighbors=15,use_rep='X_pca_harmony',random_state=0);sc.tl.umap(a,min_dist=.5,random_state=0);sc.tl.leiden(a,resolution=.5,random_state=0,key_added='subcluster_preliminary');a.uns['status']='independent subclusters; fine identities not assigned';a.write_h5ad(d/'subcluster_checkpoint.h5ad',compression='gzip')
   for c in ['donor','library','sorting','condition_Excel_original','disease status','condition_conflict']:
    if c in a.obs:pd.crosstab(a.obs.subcluster_preliminary,a.obs[c],dropna=False).to_csv(d/('subcluster_by_'+c.replace(' ','_')+'.tsv'),sep='\t')
   markers='ASGR1 CYP2C19 CYP3A4 PTH1R PTPRB STAB2 CD2 CD3D CD3E CD163 LYZ MARCO COL3A1 DCN KRT7 CFTR CD79A MS4A1 IGHM PDGFRB RBP1 LRAT ACTA2 MYH11 COL1A1 RGS5 C1QA CSF1R FCGR3A S100A8 S100A9 VCAM1 SAA1 DTNA'.split();rows=[];symbols=a.raw.var.gene_symbol_original.astype(str).to_numpy()
   for marker in markers:
    ix=np.flatnonzero(symbols==marker)
    if len(ix)!=1:continue
    val=a.raw.X[:,ix[0]].toarray().ravel()
    for cl,ci in a.obs.groupby('subcluster_preliminary',observed=True).indices.items():rows.append({'subcluster':cl,'marker':marker,'cells':len(ci),'mean_log1p_cp10k':float(val[ci].mean()),'fraction_detected':float((val[ci]>0).mean())})
   pd.DataFrame(rows).to_csv(d/'subcluster_marker_evidence.tsv',sep='\t',index=False);co=a.obs.copy();co['UMAP1']=a.obsm['X_umap'][:,0];co['UMAP2']=a.obsm['X_umap'][:,1];co.to_csv(d/'subcluster_coordinates_metadata.tsv.gz',sep='\t');result={'cells':a.n_obs,'donors':a.obs.donor.nunique(),'subclusters':a.obs.subcluster_preliminary.nunique(),'marker_rows':len(rows),'fine_identity_final':False};(d/'result.json').write_text(json.dumps(result));log(task,'complete',**result);del a
  except Exception:(d/'error.txt').write_text(traceback.format_exc());log(task,'failed',error=traceback.format_exc())
