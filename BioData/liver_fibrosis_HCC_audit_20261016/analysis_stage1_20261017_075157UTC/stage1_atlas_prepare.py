import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4])
from pathlib import Path
import json,datetime,traceback,numpy as np,pandas as pd,scanpy as sc,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).parent;OUT=O/'F2_atlas_preparation_v1';OUT.mkdir(exist_ok=False)
params={'HVG':3000,'HVG_flavor':'seurat_on_log1p_CP10k','PCs':30,'neighbors':15,'UMAP_min_dist':0.5,'leiden_resolution':0.5,'random_state':0,'scale_max_value':10,'Harmony_batch':'library','reconstructed_not_explicit_in_baseline':['PC count','neighbors','UMAP min_dist','Leiden resolution','HVG flavor','Scanpy implementation rather than original Seurat4.3'],'scope':'within-cohort preparation; no final identities, cross-study integration or disease effect claim'}
(OUT/'parameters.json').write_text(json.dumps(params,indent=2))
def log(g,status,**kw):
 x={'cohort':g,'status':status,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw};(OUT/'state.json').write_text(json.dumps(x));print(json.dumps(x),flush=True)
for g in ['GSE244832','GSE136103']:
 try:
  log(g,'loading');d=OUT/g;d.mkdir();a=sc.read_h5ad(O/g/'qc_pass_counts.h5ad');assert a.obs_names.is_unique
  sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,n_top_genes=3000,flavor='seurat')
  a.var.to_csv(d/'HVG_selection.tsv',sep='\t');a=a[:,a.var.highly_variable].copy();log(g,'PCA',cells=a.n_obs)
  sc.pp.scale(a,max_value=10);sc.tl.pca(a,n_comps=30,svd_solver='arpack',random_state=0)
  np.save(d/'PCA_unadjusted.npy',a.obsm['X_pca']);sc.pp.neighbors(a,n_neighbors=15,n_pcs=30,random_state=0);sc.tl.umap(a,min_dist=0.5,random_state=0);a.obsm['X_umap_unadjusted']=a.obsm['X_umap'].copy()
  log(g,'Harmony');sc.external.pp.harmony_integrate(a,key='library',basis='X_pca',adjusted_basis='X_pca_harmony',random_state=0,max_iter_harmony=20)
  sc.pp.neighbors(a,n_neighbors=15,use_rep='X_pca_harmony',random_state=0);sc.tl.umap(a,min_dist=0.5,random_state=0);sc.tl.leiden(a,resolution=0.5,random_state=0,key_added='cluster_preliminary')
  coords=a.obs.copy()
  for k in ['X_umap_unadjusted','X_umap']:
   for j in range(2):coords[k+'_'+str(j+1)]=a.obsm[k][:,j]
  coords.to_csv(d/'cell_coordinates_and_source_labels.tsv.gz',sep='\t')
  for label in ['library','donor','cluster_preliminary']:
   fig,axs=plt.subplots(1,2,figsize=(14,6));cats=pd.Categorical(a.obs[label]);colors=plt.get_cmap('tab20')(np.arange(len(cats.categories))%20)
   for ax,key,title in zip(axs,['X_umap_unadjusted','X_umap'],['Unadjusted','Harmony by library']):
    for j,cat in enumerate(cats.categories):
     sel=cats.codes==j;xy=a.obsm[key][sel];ax.scatter(xy[:,0],xy[:,1],s=0.3,color=colors[j],label=str(cat),rasterized=True)
    ax.set_title(title);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2')
   axs[1].legend(markerscale=7,fontsize=6,loc='upper left',bbox_to_anchor=(1,1));fig.suptitle(g+' | '+label+' | preparation, identities not finalized');fig.tight_layout();fig.savefig(d/('F2_preparation_UMAP_'+label+'.png'),dpi=180);plt.close(fig)
  a.uns['reconstruction_parameters']=json.dumps(params);a.uns['data_semantics']='scaled log-normalized HVG only; untouched QC counts remain in cohort qc_pass_counts.h5ad';a.write_h5ad(d/'atlas_HVG_preparation.h5ad',compression='gzip')
  result={'cells':a.n_obs,'HVG':a.n_vars,'clusters':int(a.obs.cluster_preliminary.nunique()),'status':'within_cohort_atlas_preparation_complete','final_cell_identity':False};(d/'result.json').write_text(json.dumps(result));log(g,'complete',**{k:v for k,v in result.items() if k!='status'});del a
 except Exception:log(g,'failed',error=traceback.format_exc())
