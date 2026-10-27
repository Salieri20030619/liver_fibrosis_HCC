from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
from scipy.stats import spearmanr
import os
os.environ['NUMBA_NUM_THREADS']='4'
import scanpy as sc
import numba
numba.set_num_threads(4)
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';I=P/'continuation_annotation_evidence_v1';O=P/'continuation_annotation_refined_v1';O.mkdir(exist_ok=True)
a=ad.read_h5ad(I/'core55_record_counts.h5ad');m=a.obs.copy();sets=json.loads((I/'marker_panels.json').read_text())['panels']
# Manual broad review uses complete81-marker evidence and whole-transcriptome cluster DE.
labels={0:'Hepatocyte',1:'Hepatocyte',2:'Hepatocyte',3:'Endothelial',4:'Endothelial',5:'Stromal',6:'Hepatocyte',7:'Myeloid',8:'Cholangiocyte',9:'Hepatocyte',10:'T_NK',11:'Endothelial',12:'Mural',13:'Hepatocyte_immune_mixed',14:'Hepatocyte',15:'Hepatocyte',16:'Hepatocyte_endothelial_mixed',17:'Cholangiocyte',18:'Myeloid'}
m['reviewed_broad_lineage']=m.seurat_clusters.astype(int).map(labels);m['cycling_state']=m.seurat_clusters.astype(int).isin([14,18]);m['annotation_basis']='81-marker panel plus full cluster DE; uncertain13/16 retained mixed pending subcluster review'
pd.DataFrame([{'cluster':k,'lineage':v,'decision':'mixed_retained' if 'mixed' in v else 'marker_supported','cycling_state':k in [14,18]} for k,v in labels.items()]).to_csv(O/'broad_annotation_decisions.tsv',sep='\t',index=False)
refnames=pd.read_csv(O/'external_reference_profile_correlations.tsv',sep='\t').reference_label.unique().tolist()
# Existing Harmony representation is reused for within-lineage neighborhood refinement.
h=pd.read_csv(P/'core55_Seurat_Harmony_repaired_v3/harmony_coordinates.tsv',sep='\t',index_col=0).loc[a.obs_names].to_numpy();m['refined_cluster']='';submarks=[];subevidence=[]
for label in sorted(m.reviewed_broad_lineage.unique()):
 ii=np.flatnonzero(m.reviewed_broad_lineage.eq(label));s=a[ii].copy();s.obsm['X_harmony']=h[ii];sc.pp.normalize_total(s,target_sum=1e4);sc.pp.log1p(s);sc.pp.neighbors(s,n_neighbors=20,use_rep='X_harmony');sc.tl.leiden(s,resolution=.4,key_added='subcluster',random_state=20261020);sc.tl.rank_genes_groups(s,'subcluster',method='wilcoxon',use_raw=False,n_genes=100)
 markers=sc.get.rank_genes_groups_df(s,group=None);markers['parent_lineage']=label;submarks.append(markers);m.loc[s.obs_names,'refined_cluster']=label+'_'+s.obs.subcluster.astype(str)
 for c in s.obs.subcluster.cat.categories:
  idx=s.obs.subcluster.eq(c).to_numpy();mean=np.asarray(s.X[idx].mean(0)).ravel();raw=a[ii[idx]].X;det=np.asarray((raw>0).mean(0)).ravel()
  for lin,genes in sets.items():
   for g in genes:
    if g in a.var_names:
     j=a.var_names.get_loc(g);subevidence.append(dict(parent_lineage=label,subcluster=c,n_cells=int(idx.sum()),lineage=lin,gene=g,mean_log1pCP10k=mean[j],pct_detected=det[j]))
 markers.to_csv(O/('checkpoint_'+label+'_DE.tsv.gz'),sep='\t',index=False);m.loc[s.obs_names].to_csv(O/('checkpoint_'+label+'_metadata.tsv.gz'),sep='\t');print('subclusters',label,len(ii),s.obs.subcluster.nunique(),flush=True)
pd.concat(submarks).to_csv(O/'subcluster_top100_DE.tsv.gz',sep='\t',index=False);pd.DataFrame(subevidence).to_csv(O/'subcluster_complete_marker_evidence.tsv.gz',sep='\t',index=False)
m.to_csv(O/'cell_metadata_reviewed.tsv.gz',sep='\t');v={'status':'reference_and_subclustering_complete_pending_review','cells':len(m),'records':55,'confirmed_crossstudy_unique_people':None,'broad_counts':m.reviewed_broad_lineage.value_counts().to_dict(),'reference_labels':refnames,'reference_is_ground_truth':False,'subclusters':m.refined_cluster.nunique(),'method':'Full81-marker evidence and whole-cluster DE; external donor-balanced profile correlations on2000 variable genes; per-lineage Leiden resolution0.4 on original Harmony30 coordinates; Wilcoxon top100 subcluster markers. Mixed groups retained; no forced assignment.','pollution':'Existing per-cell lineage codetection diagnostics, not doublet or ambient correction; empty droplets unavailable.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
