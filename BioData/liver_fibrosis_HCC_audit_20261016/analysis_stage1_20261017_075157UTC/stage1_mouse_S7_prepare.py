import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import gzip,csv,json,datetime,numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,scanpy as sc,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'GSE166504_S1_S7_preparation_v1';O.mkdir(exist_ok=False);source=P/'gap_download';m=pd.read_csv(source/'GSE166504_cell_metadata.20220204.tsv.gz',sep='\t');parts=m.FileName.str.extract(r'(?P<fraction>NPC|Hepatocyte)_(?P<timepoint>Chow|[0-9]+weeks)_Animal(?P<animal_number>[0-9]+)_Capture(?P<capture>[0-9]+)');m=pd.concat([m,parts],axis=1);m.index=m.FileName+'_'+m.CellID;m['source_unit']=m.fraction+'|'+m.timepoint+'|'+m.animal_number
sel=m.timepoint.isin(['15weeks','30weeks'])&m.CellType.isin(['Kupffer cells','Monocyte/Monocyte derived macrophage']);indices=np.flatnonzero(sel.to_numpy());rows=[];genes=[];targets={};totals=np.zeros(len(m));f=source/'GSE166504_cell_raw_counts.20220204.txt.gz'
def state(status,**kw):
 v={'status':status,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw};(O/'state.json').write_text(json.dumps(v));print(json.dumps(v),flush=True)
state('streaming_myeloid_15_30w',selected_cells=len(indices))
with gzip.open(f,'rt') as inp:
 header=next(csv.reader(inp,delimiter='\t'));cols=header[1:] if len(header)==len(m)+1 else header;assert cols==m.index.tolist()
 for line in inp:
  gene,rest=line.rstrip('\r\n').split('\t',1);gene=gene.strip('"');x=np.fromstring(rest,sep='\t');assert len(x)==len(m);totals+=x;rows.append(sp.csr_matrix(x[indices].reshape(1,-1)));genes.append(gene)
  if gene in ['Runx2','Plg','Pard3','Dtna']:targets[gene]=x.copy()
for gene,x in targets.items():m[gene+'_log1p_CP10k']=np.log1p(x/np.maximum(totals,1)*1e4)
group=['fraction','timepoint','animal_number','source_unit','CellType'];summ=m.groupby(group,observed=True).agg(cells=('FileName','size'),libraries=('FileName','nunique'),**{gene:(gene+'_log1p_CP10k','mean') for gene in targets}).reset_index();summ.to_csv(O/'source_unit_celltype_marker_summary.tsv',sep='\t',index=False)
X=sp.vstack(rows,format='csr').T.tocsr().astype('float32');del rows;a=ad.AnnData(X,obs=m.loc[sel].copy(),var=pd.DataFrame(index=pd.Index(genes)));a.var_names_make_unique();a=a[:,np.asarray(a.X.sum(axis=0)).ravel()>0].copy();a.write_h5ad(O/'myeloid_15_30w_source_counts.h5ad',compression='gzip');state('PCA_UMAP',cells=a.n_obs)
sc.pp.normalize_total(a,target_sum=1e4);sc.pp.log1p(a);sc.pp.highly_variable_genes(a,n_top_genes=3000,flavor='seurat');a=a[:,a.var.highly_variable].copy();sc.pp.scale(a,max_value=10);sc.tl.pca(a,n_comps=30,random_state=0,svd_solver='arpack');sc.pp.neighbors(a,n_neighbors=15,n_pcs=30,random_state=0);sc.tl.umap(a,min_dist=.5,random_state=0)
coords=a.obs.copy();coords['UMAP1']=a.obsm['X_umap'][:,0];coords['UMAP2']=a.obsm['X_umap'][:,1];coords.to_csv(O/'S7_myeloid_coordinates_metadata.tsv.gz',sep='\t');fig,axs=plt.subplots(1,2,figsize=(12,5));vmax=max(float(coords.Runx2_log1p_CP10k.max()),.01)
for ax,week in zip(axs,['15weeks','30weeks']):
 q=coords[coords.timepoint.eq(week)];im=ax.scatter(q.UMAP1,q.UMAP2,c=q.Runx2_log1p_CP10k,s=3,cmap='viridis',vmin=0,vmax=vmax,rasterized=True);ax.set_title('Runx2 | '+week);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');fig.colorbar(im,ax=ax)
fig.suptitle('S7 preparation | author-annotated myeloid cells | descriptive expression');fig.tight_layout();fig.savefig(O/'S7_Runx2_15_30w_preparation.png',dpi=180);plt.close(fig);a.uns['scope']='Source annotations; original processed-cell inclusion retained. Unadjusted joint15/30w embedding, not a lineage or treatment-effect claim; no automatic identity change';a.write_h5ad(O/'S7_embedding_checkpoint.h5ad',compression='gzip')
(O/'parameters_and_units.json').write_text(json.dumps({'baseline':'S7 15/30-week myeloid Runx2 expression; top3000 HVG/PCA/UMAP workflow','reconstructed':'Scanpy seurat flavor;30PC;15neighbors;min_dist0.5;seed0;no condition-confounded Harmony adjustment','QC':'uses public processed cells with author annotations; no additional cutoff imposed in this validation preparation','units':'Capture repeats combined within fraction/time/animal source units; cross-fraction animal identity not assumed; no longitudinal animal matching','genes_are_original_panel_candidates_not_confirmed_axis':list(targets),'differential_test_run':False},indent=2));state('complete',cells=a.n_obs,source_units=a.obs.source_unit.nunique(),targets=list(targets))
