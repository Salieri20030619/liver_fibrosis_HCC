from pathlib import Path
import time,json,numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
from scipy.io import mmwrite
P=Path(__file__).parent;I=P/'continuation_annotation_refined_v1';O=P/'continuation_network_communication_inputs_v1';O.mkdir(exist_ok=False)
for i in range(180):
 if (I/'validation.json').exists():break
 time.sleep(10)
else:raise RuntimeError('Refined metadata not ready')
m=pd.read_csv(I/'cell_metadata_reviewed.tsv.gz',sep='\t',index_col=0);a=ad.read_h5ad(P/'continuation_annotation_evidence_v1/core55_record_counts.h5ad');assert a.obs_names.tolist()==m.index.tolist();m['labels']=m.reviewed_broad_lineage;m.loc[m.labels.eq('Myeloid'),'labels']=m.loc[m.labels.eq('Myeloid'),'refined_cluster'];eligible=~m.reviewed_broad_lineage.str.contains('mixed');m['source_record']=m.study+'|'+m.donor.astype(str)
cc=O/'CellChat';cc.mkdir();pd.Series(a.var_names,name='gene').to_csv(cc/'gene_ids.tsv',sep='\t',index=False);aud=[]
for record in sorted(m.source_record.unique()):
 use=m.source_record.eq(record)&eligible;z=m[use].copy();counts=z.labels.value_counts();valid=counts[counts>=10].index;z=z[z.labels.isin(valid)];reason='eligible' if z.labels.nunique()>=2 else 'fewer_than2_supported_groups_with10cells'
 aud.append({'source_record':record,'study':m.loc[use,'study'].iloc[0] if use.any() else record.split('|')[0],'F':m.loc[m.source_record.eq(record),'F'].iloc[0],'eligible_cells':len(z),'groups':z.labels.nunique(),'status':reason})
 if reason!='eligible':continue
 q=cc/record.replace('|','__');q.mkdir();x=a[z.index].X.copy().astype(float);x=x.multiply((1e4/np.maximum(np.asarray(x.sum(1)).ravel(),1))[:,None]).tocsr();x.data=np.log1p(x.data);mmwrite(q/'normalized_expression.mtx',x.T);z.to_csv(q/'cell_metadata.tsv',sep='\t');print('CellChat input',record,len(z),flush=True)
pd.DataFrame(aud).to_csv(cc/'record_eligibility.tsv',sep='\t',index=False)
# Source-record/state-balanced macrophage/myeloid network. AdvancedF3/F4 discovery; allF projection.
n=O/'SCENIC';n.mkdir();my=m[m.reviewed_broad_lineage.eq('Myeloid')].copy();chosen=[]
for (record,state),g in my.groupby(['source_record','refined_cluster'],observed=True):chosen.extend(g.sample(min(150,len(g)),random_state=20261020).index)
my=my.loc[chosen].copy();my['regulon_role']=np.where(my.F.astype(float)>=3,'discovery_advanced_F3_F4','projection_F0_F2');my['subcluster_preliminary']=my.refined_cluster;x=a[my.index].X;keep=np.asarray((x>0).sum(0)).ravel()>=max(3,int(np.ceil(len(my)*.01)));genes=a.var_names[keep];frame=pd.DataFrame(x[:,keep].toarray(),index=my.index,columns=genes);frame.to_csv(n/'projection_raw_counts.tsv',sep='\t');frame.loc[my.regulon_role.eq('discovery_advanced_F3_F4')].to_csv(n/'discovery_raw_counts.tsv',sep='\t');my.to_csv(n/'cell_metadata.tsv',sep='\t');
v={'status':'prepared','all_source_records':55,'confirmed_crossstudy_independent_people':None,'CellChat_eligible_records':sum(z['status']=='eligible' for z in aud),'network_cells':len(my),'network_source_records':my.source_record.nunique(),'discovery_cells':int(my.regulon_role.eq('discovery_advanced_F3_F4').sum()),'discovery_source_records':my.loc[my.regulon_role.eq('discovery_advanced_F3_F4'),'source_record'].nunique(),'network_genes':len(genes),'network_states':my.refined_cluster.nunique(),'network_sampling':'at most150 cells per source-record/state, fixedseed20261020; all others retained in atlas. Genes detected>=1% projectioncells min3; rawcounts; motifreferences reused.','communication':'All supported cells, no cell subsampling; numeric myeloid states; >=10cells pergroup perrecord, >=2groups. Mixed13/16excluded from inference but retained in atlas.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
