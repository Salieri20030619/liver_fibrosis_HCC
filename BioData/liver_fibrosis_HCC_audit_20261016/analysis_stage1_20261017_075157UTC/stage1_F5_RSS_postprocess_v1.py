from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.special import rel_entr
from pyscenic.rss import regulon_specificity_scores
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
P=Path(__file__).parent;I=P/'F5_GSE244832_regulon_v1';O=I/'RSS_module_projection_v1';O.mkdir(exist_ok=False)
m=pd.read_csv(I/'cell_metadata.tsv',sep='\t',index_col=0,keep_default_na=False);a=pd.read_csv(I/'aucell_matrix.csv',index_col=0);assert a.index.equals(m.index)
d=m.regulon_role.eq('discovery_MASH');x=a.loc[d];lab=m.loc[d,'subcluster_preliminary'].astype(str);good=x.sum().gt(0);excluded=x.columns[~good].tolist();x=x.loc[:,good];rss=regulon_specificity_scores(x,lab)
for cl in rss.index:
 q=lab.eq(cl).to_numpy(dtype=float);q/=q.sum()
 for gene in rss.columns:
  p=x[gene].to_numpy();p=p/p.sum();mid=(p+q)/2;v=1-np.sqrt((rel_entr(p,mid).sum()+rel_entr(q,mid).sum())/2);assert np.isclose(v,rss.loc[cl,gene],atol=1e-12)
rss.to_csv(O/'discovery_all_subcluster_RSS.tsv',sep='\t');ranked=[]
for cl,row in rss.iterrows():
 z=pd.DataFrame({'regulon':row.index,'RSS':row.values}).sort_values(['RSS','regulon'],ascending=[False,True]);z['rank']=np.arange(1,len(z)+1);z['subcluster']=cl;ranked.append(z)
allr=pd.concat(ranked,ignore_index=True);allr.to_csv(O/'all_RSS_rankings.tsv',sep='\t',index=False);top=allr[allr['rank']<=10];top.to_csv(O/'top10_RSS_per_subcluster.tsv',sep='\t',index=False)
# Donor-specific RSS is a descriptive sensitivity; absent clusters are missing, never zero.
ds=[]
for donor,idx in m.loc[d].groupby('donor').groups.items():
 y=a.loc[idx,x.columns];labels=m.loc[idx,'subcluster_preliminary'].astype(str)
 if labels.nunique()<2:continue
 y=y.loc[:,y.sum().gt(0)];s=regulon_specificity_scores(y,labels);s.index.name='subcluster';z=s.stack().rename('RSS').reset_index();z.columns=['subcluster','regulon','RSS'];z['donor']=donor;ds.append(z)
don=pd.concat(ds,ignore_index=True);don.to_csv(O/'donor_subcluster_RSS.tsv',sep='\t',index=False);don.groupby(['subcluster','regulon']).RSS.agg(['mean','median','count']).to_csv(O/'equal_donor_RSS_summary.tsv',sep='\t')
fig,axes=plt.subplots(2,4,figsize=(18,9))
for ax,(cl,g) in zip(axes.flat,top.groupby('subcluster',sort=True)):
 g=g.sort_values('RSS');ax.barh(g.regulon,g.RSS,color='#347f8b');ax.set_xlim(0,1);ax.set_title('Numeric cluster '+str(cl));ax.tick_params(axis='y',labelsize=7);ax.set_xlabel('RSS (discovery cells)')
fig.suptitle('Motif-pruned regulons: top10 RSS per reviewed macrophage-like cluster');fig.tight_layout();fig.savefig(O/'top10_RSS.png',dpi=180);plt.close(fig)
mods=pd.read_csv(I/'CSI_review_v1/discovery_cells_module_assignments.tsv',sep='\t');activity=pd.DataFrame({str(k):a.loc[:,g.regulon].mean(axis=1) for k,g in mods.groupby('module_reconstructed')});activity.to_csv(O/'all_cell_module_mean_AUCell.tsv.gz',sep='\t')
# Coordinates are reused from the reviewed myeloid atlas; no graph or root is refitted.
coords=[c for c in m.columns if c.lower().replace('_','') in ['umap1','umap2']]
projection={'available':len(coords)==2,'coordinate_columns':coords,'type':'continuous activity projection, not a kernel-density estimate'}
if len(coords)==2:
 with PdfPages(O/'all_modules_existing_UMAP.pdf') as pdf:
  for module in activity.columns:
   fig,ax=plt.subplots(figsize=(7,6));z=ax.scatter(pd.to_numeric(m[coords[0]]),pd.to_numeric(m[coords[1]]),c=activity[module],s=3,cmap='viridis',rasterized=True);fig.colorbar(z,ax=ax,label='Mean AUCell across module regulons');ax.set_title('Reconstructed module '+module+' | same-cohort projection');ax.set_xlabel(coords[0]);ax.set_ylabel(coords[1]);fig.tight_layout();pdf.savefig(fig);plt.close(fig)
v={'status':'complete','RSS_source':'pySCENIC0.12.1 official regulon_specificity_scores','discovery_cells':len(x),'discovery_donors':m.loc[d,'donor'].nunique(),'regulons_scored':x.shape[1],'zero_total_AUC_regulons_excluded':excluded,'all_RSS_independently_verified_by_KL_formula':True,'disease_significance_tests':False,'UMAP_projection':projection,'independence':'same discovery cohort; not independent validation','candidate_selection':'all rankings retained; top10 display only, no forced TF'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
