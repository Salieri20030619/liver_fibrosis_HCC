from pathlib import Path
import time,json,numpy as np,pandas as pd
from scipy.stats import spearmanr
from scipy.spatial.distance import jensenshannon
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;I=P/'continuation_SCENIC_expanded_v1';M=P/'continuation_network_communication_inputs_v1/SCENIC';C=P/'continuation_SCENIC_candidates_v1';O=P/'continuation_network_source_sensitivity_v1';O.mkdir(exist_ok=False)
for i in range(1080):
 if (C/'validation.json').exists():break
 if (I/'error.txt').exists():raise RuntimeError('SCENIC upstreamfailed')
 time.sleep(10)
else:raise RuntimeError('Upstream candidatepending')
a=pd.read_csv(I/'aucell_matrix.csv',index_col=0);m=pd.read_csv(M/'cell_metadata.tsv',sep='\t',index_col=0);assert a.index.equals(m.index);dis=m.regulon_role.eq('discovery_advanced_F3_F4');balanced=a.loc[dis].groupby([m.loc[dis,'source_record'],m.loc[dis,'refined_cluster']]).mean();states=balanced.index.get_level_values(1);rss=[]
for state in sorted(set(states)):
 target=np.asarray(states==state,float);target/=target.sum();scores=[]
 for reg in balanced:
  v=balanced[reg].to_numpy();score=1-jensenshannon(v/v.sum(),target) if v.sum()>0 else np.nan;scores.append((reg,score))
 for rank,(reg,score) in enumerate(sorted(scores,key=lambda x:-x[1] if np.isfinite(x[1]) else np.inf),1):rss.append({'state':state,'regulon':reg,'record_state_balanced_RSS':score,'rank':rank})
rss=pd.DataFrame(rss);rss.to_csv(O/'record_state_balanced_RSS.tsv',sep='\t',index=False);original=pd.read_csv(C/'top10_RSS_per_state.tsv',sep='\t');overlap=[]
for state in sorted(set(states)):
 left=set(original.loc[original.subcluster.eq(state),'regulon']);right=set(rss.loc[rss.state.eq(state)&rss['rank'].le(10),'regulon']);overlap.append({'state':state,'top10_overlap':len(left&right),'original_top10_n':len(left),'balanced_top10_n':len(right)})
pd.DataFrame(overlap).to_csv(O/'RSS_sampling_sensitivity_overlap.tsv',sep='\t',index=False)
mean=a.groupby(m.source_record).mean();meta=m.groupby('source_record')[['study','F']].first();rows=[]
for study,ids in [('pooled',meta.index)]+[(name,g.index) for name,g in meta.groupby('study')]:
 for reg in mean:
  rho,p=spearmanr(mean.loc[ids,reg],meta.loc[ids,'F']);rows.append({'study':study,'regulon':reg,'TF':reg.split('(')[0],'records':len(ids),'F_rho':rho,'p_descriptive':p})
res=pd.DataFrame(rows);res['BH_q_descriptive']=np.nan
for study,g in res.groupby('study'):
 good=g.p_descriptive.notna();res.loc[g.index[good],'BH_q_descriptive']=multipletests(g.loc[good,'p_descriptive'],method='fdr_bh')[1]
res.to_csv(O/'within_source_regulon_F_associations.tsv',sep='\t',index=False);candidates=pd.read_csv(C/'F5G_candidates.tsv',sep='\t').TF.unique();res[res.TF.isin(candidates)].to_csv(O/'candidate_within_source_F_associations.tsv',sep='\t',index=False);v={'status':'completed','records':len(mean),'discovery_record_state_rows':len(balanced),'regulons':a.shape[1],'RSS_top10_overlap':overlap,'candidate_selection_changed':False,'scope':'Sensitivity only; originalcellbasedCSI/RSS candidates retained. Within-source correlations avoid attributingall source differences toF. Recordsarenotproven globallyindependent; discoveryandprojectionarenotindependentvalidation.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
