from pathlib import Path
import sys,subprocess,os,time,json
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';V=D/'F5_setup_v1/pyscenic_venv/bin/python'
if str(V)!=sys.executable:
 raise SystemExit(subprocess.call([str(V),__file__]))
import numpy as np,pandas as pd
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import pdist,jensenshannon
from scipy.stats import spearmanr
from pyscenic.cli.utils import load_signatures
I=P/'continuation_SCENIC_expanded_v1';M=P/'continuation_network_communication_inputs_v1/SCENIC';O=P/'continuation_SCENIC_candidates_v1';O.mkdir(exist_ok=False)
for i in range(720):
 if (I/'AUCell_COMPLETE.txt').exists():break
 if (I/'error.txt').exists():raise RuntimeError((I/'error.txt').read_text())
 time.sleep(10)
else:raise RuntimeError('SCENIC upstream incomplete after2hours')
a=pd.read_csv(I/'aucell_matrix.csv',index_col=0);m=pd.read_csv(M/'cell_metadata.tsv',sep='\t',index_col=0);assert a.index.equals(m.index);regs=load_signatures(str(I/'regulons.csv'));edges=[dict(regulon=r.name,TF=r.transcription_factor,target=g,weight=w) for r in regs for g,w in r.gene2weight.items()];pd.DataFrame(edges).to_csv(O/'all_regulon_targets.tsv',sep='\t',index=False)
dis=m.regulon_role.eq('discovery_advanced_F3_F4');train=a.loc[dis];good=train.std().gt(1e-14);train=train.loc[:,good];corr=np.corrcoef(train.to_numpy().T);n=len(corr);assert n>=3;csi=np.eye(n)
for i in range(n):
 for j in range(i+1,n):
  mask=np.ones(n,bool);mask[[i,j]]=False;csi[i,j]=csi[j,i]=np.mean((corr[i,mask]<corr[i,j])&(corr[j,mask]<corr[i,j]))
pd.DataFrame(csi,index=train.columns,columns=train.columns).to_csv(O/'discovery_CSI.tsv',sep='\t');Z=linkage(pdist(csi),method='complete');mods=pd.Series(fcluster(Z,t=.7*Z[:,2].max(),criterion='distance'),index=train.columns,name='module');mods.rename_axis('regulon').to_csv(O/'module_assignments.tsv',sep='\t');rss=[]
for cl in sorted(m.loc[dis,'refined_cluster'].unique()):
 target=m.loc[dis,'refined_cluster'].eq(cl).astype(float).to_numpy();target/=target.sum();scores=[]
 for reg in train:
  prob=train[reg].to_numpy();prob=prob/prob.sum();scores.append((reg,1-jensenshannon(prob,target)))
 for rank,(reg,score) in enumerate(sorted(scores,key=lambda t:t[1],reverse=True),1):rss.append({'subcluster':cl,'regulon':reg,'RSS':score,'rank':rank})
rss=pd.DataFrame(rss);rss.to_csv(O/'all_RSS.tsv',sep='\t',index=False);rss[rss['rank']<=10].to_csv(O/'top10_RSS_per_state.tsv',sep='\t',index=False)
state_mean=train.groupby(m.loc[dis,'refined_cluster']).mean();hits=[];scores=[]
for mod in sorted(mods.unique()):
 names=mods.index[mods.eq(mod)];mact=state_mean[names].mean(1);cl=mact.idxmax()
 for state,val in mact.items():scores.append({'module':int(mod),'state':state,'mean_AUCell':val,'selected_highest_state':state==cl})
 for _,row in rss[(rss.subcluster==cl)&(rss['rank']<=10)&rss.regulon.isin(names)].iterrows():hits.append({'module':int(mod),'state':cl,'regulon':row.regulon,'TF':row.regulon.split('(')[0],'RSS':row.RSS,'rank':int(row['rank'])})
pd.DataFrame(scores).to_csv(O/'module_state_scores.tsv',sep='\t',index=False);h=pd.DataFrame(hits,columns=['module','state','regulon','TF','RSS','rank']);h.to_csv(O/'F5G_candidates.tsv',sep='\t',index=False);a.groupby(m.source_record).mean().to_csv(O/'source_record_mean_AUCell.tsv',sep='\t');a.groupby([m.source_record,m.refined_cluster]).mean().to_csv(O/'source_record_state_mean_AUCell.tsv',sep='\t')
# Separate persisted-output reader for real downstream candidate association.
targets=pd.read_csv(O/'F5G_candidates.tsv',sep='\t').TF.unique();bulk=pd.read_csv(D/'GSE164760_gene_expression_median_v1.tsv.gz',sep='\t',index_col=0);fib=pd.read_csv(D/'F6E_GSE164760_fibrosis_score_v1/official_ssGSEA_scores.tsv',sep='\t',index_col=0);assert bulk.columns.equals(fib.columns);rows=[]
for gene in targets:
 if gene in bulk.index:
  rho,p=spearmanr(bulk.loc[gene],fib.iloc[0]);rows.append({'TF':gene,'n':len(bulk.columns),'rho':rho,'p_descriptive':p})
pd.DataFrame(rows).to_csv(O/'downstream_F6D_bulk_fibrosis_correlations.tsv',sep='\t',index=False);v={'status':'completed','regulons':len(regs),'target_edges':len(edges),'projection_cells':len(a),'source_records':m.source_record.nunique(),'discovery_cells':int(dis.sum()),'discovery_records':m.loc[dis,'source_record'].nunique(),'states':m.refined_cluster.nunique(),'candidate_TFs':len(targets),'downstream_matched_TFs':len(rows),'selection':'EachCSI module assigned highestcellmeanAUCell state, intersect thatstate top10RSS; modulecut70%maxheight reconstructed; noDTNAoraxisassumption','limits':'F3/F4 crosssectional discovery; source-record sampling; globallyindependentpersonsnotconfirmed; pooledG164760correlation descriptive only, notincidentHCC.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
