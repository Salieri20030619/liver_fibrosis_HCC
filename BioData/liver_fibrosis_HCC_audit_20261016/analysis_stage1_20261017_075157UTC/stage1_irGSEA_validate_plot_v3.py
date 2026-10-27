from pathlib import Path
import json,numpy as np,pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;source=P/'F4_irGSEA_official_v2';O=P/'F4_irGSEA_validated_v3';O.mkdir(exist_ok=False);allresults=[]
for g in ['GSE136103','GSE244832']:
 d=source/g;meta=pd.read_csv(P/'F4_AUCell_official_v1'/g/'cell_metadata.tsv',sep='\t',index_col=0);out=O/g;out.mkdir()
 for sort,z in meta.groupby('sorting'):
  dm=z[['donor','group']].drop_duplicates().set_index('donor');assert dm.index.is_unique;target='fibrotic' if g=='GSE136103' else 'F3_F4';reference='healthy' if g=='GSE136103' else 'F0';rows=[];mat=[]
  for method in ['AUCell','UCell','singscore','ssgsea']:
   f=d/(str(sort)+'_'+method+'_donor_scores.tsv');scores=pd.read_csv(f,sep='\t',index_col=0);assert scores.shape[0]==50 and set(scores.columns)==set(dm.index) and np.isfinite(scores).all().all();records=[]
   for gene,v in scores.iterrows():
    a=v[dm.index[dm.group.eq(target)]];b=v[dm.index[dm.group.eq(reference)]];records.append({'pathway':gene,'method':method,'target':target,'reference':reference,'n_target':len(a),'n_reference':len(b),'mean_difference':a.mean()-b.mean(),'p':mannwhitneyu(a,b,alternative='two-sided',method='asymptotic').pvalue})
   t=pd.DataFrame(records);p=t.p.to_numpy();order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);t['BH_q']=q;rows.append(t);mat.append(pd.Series((np.sign(t.mean_difference)*-np.log10(t.p.clip(lower=1e-300))).to_numpy(),index=t.pathway.to_numpy(),name=method))
  t=pd.concat(rows);key=str(sort).replace('+','positive').replace('-','negative');t.to_csv(out/(key+'_donor_tests.tsv'),sep='\t',index=False);m=pd.concat(mat,axis=1);fig,ax=plt.subplots(figsize=(9,13));v=max(abs(m.to_numpy()).max(),1);im=ax.imshow(m,cmap='RdBu_r',aspect='auto',vmin=-v,vmax=v);ax.set_yticks(range(len(m)),m.index.str.replace('HALLMARK-',''),fontsize=7);ax.set_xticks(range(4),m.columns,rotation=25);ax.set_title(g+' '+str(sort)+' | '+target+' vs '+reference+'\nFour official methods; donor-level companion tests');fig.colorbar(im,ax=ax,label='Signed -log10(raw p); not RRA significance');fig.tight_layout();fig.savefig(out/(key+'_four_method_donor_heatmap.png'),dpi=160);plt.close(fig)
  integrations=[str(f.name) for f in d.glob(str(sort)+'*_integration.tsv')];errors=[f.read_text() for f in d.glob(str(sort)+'*_integration_error.txt')];allresults.append({'cohort':g,'sorting':str(sort),'donors':len(dm),'pathways_per_method':50,'BH_q_below05_by_method':t.groupby('method').BH_q.apply(lambda x:int((x<.05).sum())).to_dict(),'official_integration_files':integrations,'official_integration_errors':errors})
(O/'validation.json').write_text(json.dumps(allresults,indent=2));print(json.dumps(allresults))
