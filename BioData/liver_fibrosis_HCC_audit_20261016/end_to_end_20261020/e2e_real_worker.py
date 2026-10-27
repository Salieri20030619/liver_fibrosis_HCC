from pathlib import Path
import json,hashlib,traceback,time
import numpy as np,pandas as pd
from scipy.stats import spearmanr
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';I=D/'F5_GSE244832_regulon_v1'
O=P/'F5_to_F6_actual';O.mkdir(exist_ok=False)
try:
 a=pd.read_csv(I/'aucell_matrix.csv',index_col=0);m=pd.read_csv(I/'cell_metadata.tsv',sep='\t',index_col=0)
 assert a.index.equals(m.index)
 rss=pd.read_csv(I/'RSS_module_projection_v1/top10_RSS_per_subcluster.tsv',sep='\t')
 modules=pd.read_csv(I/'CSI_review_v1/discovery_cells_module_assignments.tsv',sep='\t')
 d=m.regulon_role.eq('discovery_MASH'); assert d.sum()==2022
 # Equal donor means prevent large cell libraries from driving candidate module ranking.
 z=a.loc[d].groupby([m.loc[d,'donor'],m.loc[d,'subcluster_preliminary']]).mean().groupby(level=1).mean()
 scores=[];hits=[]
 for cl,row in z.iterrows():
  values={int(mod):float(row[g.regulon].mean()) for mod,g in modules.groupby('module_reconstructed')}
  best=max(values,key=values.get)
  for mod,v in values.items():scores.append({'subcluster':str(cl),'module':mod,'mean_AUC_equal_donor':v,'selected_highest_mean':mod==best})
  rs=rss[rss.subcluster.astype(str).eq(str(cl))]
  for _,v in rs.merge(modules,on='regulon').query('module_reconstructed == @best').iterrows():hits.append({'subcluster':str(cl),'module':best,'regulon':v.regulon,'TF':v.regulon.split('(')[0],'RSS':v.RSS,'RSS_rank':v['rank']})
 pd.DataFrame(scores).to_csv(O/'module_scores.tsv',sep='\t',index=False)
 h=pd.DataFrame(hits);assert len(h)>0;h.to_csv(O/'F5G_candidates.tsv',sep='\t',index=False)
 # A separate downstream reader consumes the persisted upstream candidate file.
 targets=pd.read_csv(O/'F5G_candidates.tsv',sep='\t').TF.unique()
 bulk=pd.read_csv(D/'GSE164760_gene_expression_median_v1.tsv.gz',sep='\t',index_col=0)
 fibrosis=pd.read_csv(D/'F6E_GSE164760_fibrosis_score_v1/official_ssGSEA_scores.tsv',sep='\t',index_col=0)
 assert bulk.columns.equals(fibrosis.columns)
 results=[]
 for tf in targets:
  if tf not in bulk.index:continue
  rho,p=spearmanr(bulk.loc[tf],fibrosis.iloc[0]);results.append({'TF':tf,'n':len(bulk.columns),'rho':rho,'p_descriptive_only':p})
 pd.DataFrame(results).to_csv(O/'F6D_candidate_fibrosis_correlations.tsv',sep='\t',index=False)
 v={'status':'completed','discovery_cells':2022,'discovery_donors':int(m.loc[d,'donor'].nunique()),'projection_cells':len(m),'candidate_rows':len(h),'candidate_TFs':len(targets),'downstream_TFs_matched':len(results),'bulk_samples':bulk.shape[1],'adaptation':'Original module/RSS intersection; module selection reconstructed as highest equal-donor mean AUC for every numeric subcluster. No DTNA or RUNX2 prespecified. Correlation is descriptive pooled cross-sectional, not incident HCC or causal evidence.','input_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [I/'RSS_module_projection_v1/top10_RSS_per_subcluster.tsv',I/'CSI_review_v1/discovery_cells_module_assignments.tsv',O/'F5G_candidates.tsv']}}
 (O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
