from pathlib import Path
import time,json,traceback
import pandas as pd,numpy as np,anndata as ad
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;O=P/'new_atlas_to_TF_validation_v2';O.mkdir(exist_ok=False);atlas=P/'core55_marker_DE_annotation_v5'
try:
 for i in range(80):
  if (atlas/'cell_metadata_with_annotation.tsv').exists():break
  if (atlas/'error.txt').exists():raise RuntimeError('Atlas failed; metadata unavailable, see upstream error')
  print('WAIT real atlas annotations',i,flush=True);time.sleep(30)
 else:raise RuntimeError('Atlas annotations not ready; no proxy labels used')
 m=pd.read_csv(atlas/'cell_metadata_with_annotation.tsv',sep='\t',index_col=0);a=ad.read_h5ad(P/'GSE202379_core39_QC/core39_counts.h5ad');a.obs_names='GSE202379|'+a.obs_names
 m=m.loc[a.obs_names];assert len(m)==a.n_obs and m.donor.nunique()==39
 tf=pd.read_csv(P/'F5_to_F6_original_direction_v2/F5G_candidates.tsv',sep='\t').TF.unique();features=[g for g in tf if g in a.var_names];rows=[]
 for (donor,label),ids in m.groupby(['donor','provisional_celltype']).groups.items():
  if label not in ['Myeloid','Mesenchymal_stromal','Mesenchymal_mural']:continue
  ix=a.obs_names.get_indexer(ids);raw=np.asarray(a.X[ix].sum(axis=0)).ravel();total=raw.sum()
  for g in features:rows.append({'donor':donor,'F':int(m.loc[ids,'F'].iloc[0]),'celltype':label,'cells':len(ids),'gene':g,'pseudobulk_log2CPM':float(np.log2(raw[a.var_names.get_loc(g)]/total*1e6+1))})
 z=pd.DataFrame(rows);assert len(z)>0;z.to_csv(O/'donor_celltype_TF_expression.tsv',sep='\t',index=False);stats=[]
 for (ct,g),v in z[z.cells>=30].groupby(['celltype','gene']):
  if len(v)<5:continue
  rho,p=spearmanr(v.F,v.pseudobulk_log2CPM);stats.append({'celltype':ct,'gene':g,'donors':len(v),'rho_F':rho,'p':p})
 q=pd.DataFrame(stats);q['BH_q']=np.nan
 for ct,ids in q.groupby('celltype').groups.items():
  good=q.loc[ids,'p'].notna();use=q.loc[ids].index[good];q.loc[use,'BH_q']=multipletests(q.loc[use,'p'],method='fdr_bh')[1]
 q.to_csv(O/'source_heldout_F_correlations.tsv',sep='\t',index=False)
 v={'status':'completed','upstream':'core55 Seurat/Harmony clusters with marker-DE-supported v5 sidecar annotation physically reread','cohort':'GSE20237939core directF donors; candidateTFs frozen fromG2448328MASH discoverydonors','TFs_matched':len(features),'tests':len(q),'q_below05':int(q.BH_q.lt(.05).sum()),'minimum_cells_per_donor_celltype':30,'limitations':'Provisional coarse annotations; expression association notregulonactivity validation. Not prospective clinical validation; no reselection of bulkmodel features.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
