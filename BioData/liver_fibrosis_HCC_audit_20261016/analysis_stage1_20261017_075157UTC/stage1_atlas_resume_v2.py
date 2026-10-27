import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4])
from pathlib import Path
import json,datetime,traceback,numpy as np,pandas as pd,anndata as ad,scanpy as sc,harmonypy as hm,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'F2_atlas_preparation_v2';O.mkdir(exist_ok=False)
def log(g,status,**kw):
 x={'cohort':g,'status':status,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw};(O/'state.json').write_text(json.dumps(x));print(json.dumps(x),flush=True)
def harmony(pc,obs,**kw):
 z=np.asarray(hm.run_harmony(pc,obs,'library',random_state=0,**kw).Z_corr)
 if z.shape==pc.shape:return z
 if z.T.shape==pc.shape:return z.T
 raise ValueError(('unexpected Harmony dimensions',z.shape,pc.shape))
rng=np.random.default_rng(0);tiny=rng.normal(size=(120,8));test=harmony(tiny,pd.DataFrame({'library':['a']*60+['b']*60}),nclust=4,max_iter_harmony=2)
assert test.shape==tiny.shape and np.isfinite(test).all();(O/'compatibility_smoke_test.json').write_text(json.dumps({'passed':True,'shape':list(test.shape),'fix':'Check output dimensions; avoid unconditional transpose in incompatible Scanpy wrapper'}))
for g in ['GSE244832','GSE136103']:
 try:
  d=O/g;d.mkdir();src=ad.read_h5ad(P/g/'qc_pass_counts.h5ad',backed='r');obs=src.obs.copy();src.file.close();pc=np.load(P/'F2_atlas_preparation_v1'/g/'PCA_unadjusted.npy');assert pc.shape==(len(obs),30)
  a=ad.AnnData(obs=obs);a.obsm['X_pca']=pc;log(g,'resume_saved_PCA',cells=len(obs));z=harmony(pc,obs,max_iter_harmony=20);np.save(d/'PCA_harmony.npy',z);a.obsm['X_pca_harmony']=z
  for basis,name in [('X_pca','X_umap_unadjusted'),('X_pca_harmony','X_umap_harmony')]:
   log(g,'neighbors_UMAP',basis=basis);sc.pp.neighbors(a,n_neighbors=15,use_rep=basis,random_state=0);sc.tl.umap(a,min_dist=.5,random_state=0);a.obsm[name]=a.obsm['X_umap'].copy()
  sc.tl.leiden(a,resolution=.5,random_state=0,key_added='cluster_preliminary');coords=a.obs.copy()
  for name in ['X_umap_unadjusted','X_umap_harmony']:
   coords[name+'_1']=a.obsm[name][:,0];coords[name+'_2']=a.obsm[name][:,1]
  coords.to_csv(d/'cell_coordinates_and_source_labels.tsv.gz',sep='\t');a.uns['status']='embedding preparation only; no finalized identities; original counts remain separate';a.uns['parameters']='Reuse original 3000HVG 30PC checkpoint; neighbors15, UMAPmin_dist0.5, Leiden0.5, random_state0; explicit reconstructed parameters';a.write_h5ad(d/'embedding_checkpoint.h5ad',compression='gzip')
  for c in ['library','donor','sorting','condition_metadata_original','condition_Excel_original','disease status','disease_comparison_eligible']:
   if c in a.obs:pd.crosstab(a.obs.cluster_preliminary,a.obs[c],dropna=False).to_csv(d/('cluster_by_'+c.replace(' ','_')+'.tsv'),sep='\t')
  for label in ['library','donor','cluster_preliminary']:
   fig,axs=plt.subplots(1,2,figsize=(14,6));cats=pd.Categorical(a.obs[label]);colors=plt.get_cmap('tab20')(np.arange(len(cats.categories))%20)
   for ax,key in zip(axs,['X_umap_unadjusted','X_umap_harmony']):
    for j,cat in enumerate(cats.categories):
     xy=a.obsm[key][cats.codes==j];ax.scatter(xy[:,0],xy[:,1],s=.3,color=colors[j],label=str(cat),rasterized=True)
    ax.set_title(key)
   axs[1].legend(markerscale=7,fontsize=6,loc='upper left',bbox_to_anchor=(1,1));fig.suptitle(g+' | '+label+' | preliminary');fig.tight_layout();fig.savefig(d/('F2_preparation_'+label+'.png'),dpi=160);plt.close(fig)
  log(g,'marker_evidence');src=ad.read_h5ad(P/g/'qc_pass_counts.h5ad',backed='r');markers=['ALB','ASGR1','CYP3A4','PECAM1','VWF','STAB2','DCN','COL3A1','COL1A1','CD3D','CD3E','CD2','MS4A1','CD79A','CD79B','LYZ','CD163','MARCO','KRT7','KRT19','EPCAM','NKG7','GNLY','KLRD1','RGS5','MYH11','ACTA2'];syms=src.var.gene_symbol_original.astype(str).to_numpy();indices=[int(np.flatnonzero(syms==m)[0]) for m in markers if (syms==m).sum()==1];valid=[syms[i] for i in indices];values=np.empty((src.n_obs,len(indices)),dtype='float32');detected=np.empty(values.shape,dtype=bool)
  for start in range(0,src.n_obs,5000):
   end=min(start+5000,src.n_obs);xx=src.X[start:end,:];total=np.asarray(xx.sum(axis=1)).ravel();mm=xx[:,indices].toarray();values[start:end]=np.log1p(mm/np.maximum(total[:,None],1)*1e4);detected[start:end]=mm>0
  src.file.close();rows=[]
  for cl,ix in a.obs.groupby('cluster_preliminary',observed=True).indices.items():
   for j,m in enumerate(valid):rows.append({'cluster':cl,'marker':m,'cells':len(ix),'mean_log1p_cp10k':float(values[ix,j].mean()),'fraction_detected':float(detected[ix,j].mean())})
  pd.DataFrame(rows).to_csv(d/'cluster_marker_evidence.tsv',sep='\t',index=False);result={'cells':a.n_obs,'clusters':a.obs.cluster_preliminary.nunique(),'marker_rows':len(rows),'final_cell_identity':False,'status':'preparation_and_review_evidence_complete'};(d/'result.json').write_text(json.dumps(result));log(g,'complete',cells=a.n_obs);del a,values,detected
 except Exception:log(g,'failed',error=traceback.format_exc())
