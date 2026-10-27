from pathlib import Path
import json,hashlib,itertools
import numpy as np,pandas as pd
from scipy.cluster.hierarchy import linkage,leaves_list,fcluster,dendrogram
from scipy.spatial.distance import pdist
from pyscenic.cli.utils import load_signatures
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;I=P/'F5_GSE244832_regulon_v1';O=I/'CSI_review_v1';O.mkdir(exist_ok=False)
m=pd.read_csv(I/'cell_metadata.tsv',sep='\t',index_col=0,keep_default_na=False);a=pd.read_csv(I/'aucell_matrix.csv',index_col=0);assert a.index.tolist()==m.index.tolist() and a.columns.is_unique and np.isfinite(a).all().all() and (a>=0).all().all()
regs=load_signatures(str(I/'regulons.csv'));assert set(a.columns)=={r.name for r in regs};edges=[]
for r in regs:
 for gene,weight in r.gene2weight.items():edges.append({'regulon':r.name,'TF':r.transcription_factor,'target':gene,'weight':weight})
pd.DataFrame(edges).to_csv(O/'all_motif_pruned_regulon_targets.tsv',sep='\t',index=False);(O/'all_regulons.gmt').write_text('\n'.join(r.name+'\tpySCENIC_motif_pruned\t'+'\t'.join(r.genes) for r in regs)+'\n')
discovery=m.regulon_role.eq('discovery_MASH');train=a.loc[discovery];donor=a.groupby(m.donor,sort=True).mean();donor.to_csv(O/'donor_mean_AUCell.tsv',sep='\t');meta=m.groupby('donor')[['condition_Excel_original','F_stage_from_Figure4C','regulon_role']].first();meta.to_csv(O/'donor_metadata.tsv',sep='\t');sub=a.groupby([m.donor,m.subcluster_preliminary],sort=True).mean();sub.to_csv(O/'donor_subcluster_mean_AUCell.tsv',sep='\t');checks=[]
def CSI(frame,name):
 good=frame.std(axis=0).gt(1e-14);excluded=frame.columns[~good].tolist();z=frame.loc[:,good];assert z.shape[1]>=3;corr=np.corrcoef(z.to_numpy().T);assert np.isfinite(corr).all();n=len(corr);csi=np.eye(n)
 for i in range(n):
  for j in range(i+1,n):
   mask=np.ones(n,dtype=bool);mask[[i,j]]=False;v=np.mean((corr[i,mask]<corr[i,j])&(corr[j,mask]<corr[i,j]));csi[i,j]=csi[j,i]=v
 # Independent explicit enumeration of all pairs for this small regulon matrix.
 for i in range(n):
  for j in range(i+1,n):assert np.isclose(csi[i,j],sum(corr[i,k]<corr[i,j] and corr[j,k]<corr[i,j] for k in range(n) if k not in [i,j])/(n-2))
 assert np.allclose(csi,csi.T) and ((csi>=0)&(csi<=1)).all();pd.DataFrame(corr,index=z.columns,columns=z.columns).to_csv(O/(name+'_Pearson.tsv'),sep='\t');pd.DataFrame(csi,index=z.columns,columns=z.columns).to_csv(O/(name+'_CSI.tsv'),sep='\t');Z=linkage(pdist(csi,metric='euclidean'),method='complete');np.save(O/(name+'_linkage.npy'),Z);order=leaves_list(Z);modules=fcluster(Z,t=.7*Z[:,2].max(),criterion='distance');pd.DataFrame({'regulon':z.columns,'module_reconstructed':modules}).to_csv(O/(name+'_module_assignments.tsv'),sep='\t',index=False);es=[{'regulon1':z.columns[i],'regulon2':z.columns[j],'CSI':csi[i,j]} for i,j in itertools.combinations(range(n),2) if csi[i,j]>.7];pd.DataFrame(es,columns=['regulon1','regulon2','CSI']).to_csv(O/(name+'_CSI_above_0_7_edges.tsv'),sep='\t',index=False)
 fig,ax=plt.subplots(figsize=(9,8));im=ax.imshow(csi[order][:,order],vmin=0,vmax=1,cmap='viridis');ax.set_title(name+' | all nonconstant regulons\nCSI; complete-linkage Euclidean ordering');ax.set_xlabel('Regulon (hierarchical order)');ax.set_ylabel('Regulon (hierarchical order)');fig.colorbar(im,ax=ax,label='CSI');fig.tight_layout();fig.savefig(O/(name+'_CSI_heatmap.png'),dpi=180);plt.close(fig)
 if name=='discovery_cells':
  fig,ax=plt.subplots(figsize=(14,5));dendrogram(Z,labels=z.columns.to_numpy(),leaf_rotation=90,leaf_font_size=5,ax=ax);ax.axhline(.7*Z[:,2].max(),color='red',linestyle='--',label='Reconstructed module cut:70%max height');ax.set_ylabel('Euclidean complete-linkage distance');ax.legend();fig.tight_layout();fig.savefig(O/'regulon_hierarchical_tree.png',dpi=180,bbox_inches='tight');plt.close(fig)
 checks.append({'matrix':name,'observations':len(frame),'nonconstant_regulons':n,'excluded_constant_regulons':excluded,'modules_reconstructed':len(set(modules)),'CSI_edges_above_0_7':len(es),'all_pair_CSI_independently_recomputed':True});return pd.Series(modules,index=z.columns),order
mods,order=CSI(train,'discovery_cells');CSI(donor.loc[meta.regulon_role.eq('discovery_MASH')],'discovery_donor_means')
module_activity=pd.DataFrame({str(k):donor.loc[:,mods.index[mods.eq(k)]].mean(axis=1) for k in sorted(mods.unique())});module_activity.to_csv(O/'donor_module_mean_AUCell.tsv',sep='\t');plot=module_activity.T;std=plot.std(axis=1).replace(0,np.nan);zs=plot.sub(plot.mean(axis=1),axis=0).div(std,axis=0).fillna(0);fig,ax=plt.subplots(figsize=(12,max(4,len(zs)*.35)));im=ax.imshow(zs,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2);ax.set_yticks(range(len(zs)),['Module '+x for x in zs.index]);ax.set_xticks(range(len(zs.columns)),[d+' '+str(meta.loc[d,'condition_Excel_original'])+' F'+str(meta.loc[d,'F_stage_from_Figure4C']) for d in zs.columns],rotation=90,fontsize=8);ax.set_title('Donor-average module activity | same-cohort projections, not external validation');fig.colorbar(im,ax=ax,label='Row z-score');fig.tight_layout();fig.savefig(O/'donor_module_activity.png',dpi=180,bbox_inches='tight');plt.close(fig)
cluster=sub.groupby(level=1).mean();cluster.to_csv(O/'equal_donor_weight_subcluster_AUCell.tsv',sep='\t');show=train.var().sort_values(ascending=False).head(50).index;z=cluster[show].T;z=z.sub(z.mean(axis=1),axis=0).div(z.std(axis=1).replace(0,np.nan),axis=0).fillna(0);fig,ax=plt.subplots(figsize=(8,max(6,len(z)*.18)));im=ax.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2);ax.set_yticks(range(len(z)),z.index,fontsize=6);ax.set_xticks(range(len(z.columns)),z.columns);ax.set_xlabel('Reviewed numeric macrophage-like cluster; no forced subtype');ax.set_title('Top50 variable discovery regulons for display only\nFull regulon matrix retained; donor-balanced cluster averages');fig.colorbar(im,ax=ax,label='Row z-score');fig.tight_layout();fig.savefig(O/'subcluster_regulon_activity.png',dpi=180,bbox_inches='tight');plt.close(fig)
v={'status':'complete','AUCell_shape':list(a.shape),'cell_ID_order_exact':True,'regulons':len(regs),'regulon_target_edges':len(edges),'CSI_checks':checks,'discovery_donors':8,'projection_donors':m.donor.nunique(),'no_disease_p_values':True,'module_cut_is_reconstruction':True,'no_original_axis_assumed':True,'regulons_csv_sha256':hashlib.sha256((I/'regulons.csv').read_bytes()).hexdigest()};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
