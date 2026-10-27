import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,scipy.sparse as sp,scipy.io,json,subprocess,traceback
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'F4_AUCell_postprocess_v2';O.mkdir(exist_ok=False);results=[]
for g in ['GSE244832','GSE136103']:
 try:
  d=O/g;d.mkdir();source=P/'F4_AUCell_official_v1'/g;obs=pd.read_csv(source/'cell_metadata.tsv',sep='\t',index_col=0)
  scores=pd.read_csv(source/'official_AUCell_scores.tsv',sep='\t',index_col=0);assert scores.index.tolist()==obs.index.tolist() and np.isfinite(scores).all().all();full=obs[['donor','sorting','group']].join(scores);donor=full.groupby(['donor','sorting','group'],observed=True)[list(scores.columns)].mean().reset_index();donor.to_csv(d/'donor_AUCell_scores.tsv',sep='\t',index=False);u=pd.read_csv(P/'F4_UCell_Hallmark_v1'/(g+'_donor_Hallmark_UCell.tsv'),sep='\t');joined=donor.merge(u,on=['donor','sorting','group'],suffixes=('_AUC','_UCell'),validate='one_to_one');cor=[{'pathway':name,'Spearman_AUCell_vs_UCell':joined[name+'_AUC'].corr(joined[name+'_UCell'],method='spearman')} for name in scores.columns];pd.DataFrame(cor).to_csv(d/'AUCell_UCell_donor_concordance.tsv',sep='\t',index=False)
  for sort,z in donor.groupby('sorting',observed=True):
   groups=sorted(z.group.unique());rows=[]
   for name in scores.columns:
    xx=z[z.group.eq(groups[0])][name];yy=z[z.group.eq(groups[1])][name];rows.append({'pathway':name,'group1':groups[0],'group2':groups[1],'n1':len(xx),'n2':len(yy),'mean_difference':xx.mean()-yy.mean(),'p_value':mannwhitneyu(xx,yy,alternative='two-sided',method='asymptotic').pvalue})
   t=pd.DataFrame(rows);p=t.p_value.to_numpy();order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);t['BH_q']=q;key=str(sort).replace('+','positive').replace('-','negative');t.to_csv(d/(key+'_pathway_tests.tsv'),sep='\t',index=False);show=['HALLMARK_GLYCOLYSIS','HALLMARK_HYPOXIA','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_COAGULATION'];fig,axs=plt.subplots(1,4,figsize=(14,4))
   for ax,name in zip(axs,show):
    vals=[z[z.group.eq(h)][name] for h in groups];ax.boxplot(vals,tick_labels=groups)
    for i,v in enumerate(vals):ax.scatter(np.full(len(v),i+1),v,s=16)
    ax.set_title(name.replace('HALLMARK_',''),fontsize=8);ax.set_ylabel('Mean official AUCell per donor')
   fig.suptitle(g+' '+str(sort)+' | AUCell companion, not irGSEA integration');fig.tight_layout();fig.savefig(d/(key+'_F4E_AUCell_companion.png'),dpi=160);plt.close(fig);results.append({'cohort':g,'sorting':str(sort),'BH_q_below05':int((q<.05).sum()),'donors':len(z)})
  (O/'state.json').write_text(json.dumps(results));print(json.dumps(results),flush=True)
 except Exception:(d/'error.txt').write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
(O/'scope.json').write_text(json.dumps({'method':'Official AUCell1.32.0; seed0; aucMaxRank ceiling5% genes; normAUC TRUE. Default rank cutoff is algorithm parameter, no new cell filter.','panels':'F4E companion; full irGSEA F4G remains separate','input':'Same eligible myeloid cells, duplicate symbols summed before ranks; public Hallmark2025.1.Hs, not original unavailable custom sets','statistical_unit':'donor x sorting; BH across50 pathways per comparison; disease labels and D9/D10 exclusions retained'},indent=2))
