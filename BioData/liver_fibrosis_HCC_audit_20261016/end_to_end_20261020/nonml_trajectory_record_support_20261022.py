from pathlib import Path
import sys,time,json
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'))
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata
O=P/'nonml_trajectory_record_support_20261022';O.mkdir(exist_ok=False);start=time.time();checks=[]
for branch,folder in [('F3','continuation_F3_CytoTRACE_serial_resume_v2'),('F5','continuation_expanded_trajectories_v1/F5')]:
 q=P/folder;i=P/'continuation_followup_inputs_v1'/branch
 m=pd.read_csv(q/'ordered_cell_metadata.tsv',sep='\t',index_col=0);de=pd.read_csv(q/'pseudotime_dynamic_genes.tsv',sep='\t',index_col=0)
 cells=(i/'cells.txt').read_text().splitlines();genes=(i/'genes.txt').read_text().splitlines();assert cells==m.index.tolist()
 x=mmread(i/'raw_counts.mtx').tocsr();total=np.asarray(x.sum(0)).ravel();ok=de.status.eq('OK');ids=de.index[ok];take=pd.Index(genes).get_indexer(ids);assert (take>=0).all()
 y=x[take].toarray().astype(float);y=np.log1p(y/np.maximum(total,1)*1e4);del x
 rows=[]
 for record,g in m.groupby('donor'):
  idx=m.index.get_indexer(g.index);t=g.Pseudotime.to_numpy();n=len(idx)
  rho=np.full(len(ids),np.nan);reason='fewer_than3_cells' if n<3 else 'constant_pseudotime' if len(np.unique(t))<2 else 'computed_where_expression_varies'
  if n>=3 and len(np.unique(t))>=2:
   z=rankdata(y[:,idx],axis=1);z-=z.mean(1,keepdims=True);tr=rankdata(t);tr-=tr.mean();den=np.sqrt((z*z).sum(1)*(tr*tr).sum());valid=den>0;rho[valid]=(z[valid]@tr)/den[valid]
  rows.append(pd.DataFrame({'gene':ids,'record':record,'study':g.study.iloc[0],'F':g.F.iloc[0],'cells':n,'rho':rho,'reason':reason}))
  print(branch,record,n,flush=True)
 rd=pd.concat(rows,ignore_index=True);rd.to_csv(O/(branch+'_record_gene_correlations_local.tsv.gz'),sep='\t',index=False)
 ag=rd.groupby(['gene','study']).agg(records=('record','size'),estimable_records=('rho','count'),median_rho=('rho','median'),positive_records=('rho',lambda v:int((v>0).sum())),negative_records=('rho',lambda v:int((v<0).sum())),zero_records=('rho',lambda v:int((v==0).sum()))).reset_index();ag['unestimable_records']=ag.records-ag.estimable_records
 ag=ag.merge(de[['qval','pval','status']],left_on='gene',right_index=True,validate='many_to_one');ag.to_csv(O/(branch+'_gene_source_support.tsv'),sep='\t',index=False)
 wide=ag.pivot(index='gene',columns='study',values='median_rho');wide['same_median_sign']=np.where(wide.notna().all(axis=1),np.sign(wide.iloc[:,0])==np.sign(wide.iloc[:,1]),np.nan);wide=wide.join(de[['qval','status']]);wide.to_csv(O/(branch+'_cross_source_support.tsv'),sep='\t')
 checks.append({'branch':branch,'cells':len(m),'records':m.donor.nunique(),'successful_dynamic_genes':len(ids),'failed_genes_excluded':int((~ok).sum()),'cell_dynamic_q005':int((de.loc[ids,'qval']<.05).sum()),'record_gene_correlations':len(rd),'finite_correlations':int(rd.rho.notna().sum()),'cross_source_median_same_sign':int(wide.same_median_sign.eq(1).sum()),'dynamic_q005_same_source_median_sign':int((wide.qval.lt(.05)&wide.same_median_sign.eq(1)).sum())})
 del y,rd
(O/'validation.json').write_text(json.dumps({'branches':checks,'seconds':time.time()-start,'method':'Within-record Spearman log1pCP10K vs frozen pooled pseudotime; no new p values; constants and <3cells remain NA; all original successful tested genes included; descriptive support, not independent donor validation or newly inferred trajectory'},indent=2));print(json.dumps(checks),flush=True)
