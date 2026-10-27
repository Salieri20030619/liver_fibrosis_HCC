from pathlib import Path
import anndata as ad,pandas as pd,numpy as np,scipy.sparse as sp,json,traceback
from scipy.io import mmwrite
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'core55_QC_counts';O.mkdir(exist_ok=False)
try:
 a=ad.read_h5ad(P/'GSE202379_core39_QC/core39_counts.h5ad');b=ad.read_h5ad(D/'GSE244832/qc_pass_counts.h5ad')
 b=b[b.obs.donor.astype(str).isin(['D'+str(i) for i in list(range(1,9))+list(range(11,19))])].copy();assert b.obs.donor.nunique()==16
 sy=b.var.gene_symbol_original.astype(str);dup=sy.duplicated(keep=False);b=b[:,~dup].copy();b.var_names=b.var.gene_symbol_original.astype(str)
 common=a.var_names.intersection(b.var_names,sort=False);assert len(common)>15000
 ma=pd.read_csv(P/'GSE202379_core39_QC_MT_NA_v2/cell_metadata.tsv',sep='\t',index_col=0);ma=ma.loc[a.obs_names];mb=pd.DataFrame({'donor':b.obs.donor.astype(str),'GSM':b.obs.GSM.astype(str),'F':pd.to_numeric(b.obs.F_stage_from_Figure4C),'study':'GSE244832','UMI':b.obs.total_counts_recomputed,'genes':b.obs.n_genes_recomputed,'percent_mt':b.obs.pct_mt_recomputed,'mt_QC_applied':True},index=b.obs_names)
 for m in [ma,mb]:m['donor']=m.study+'|'+m.donor;m.index=m.study+'|'+m.index.astype(str)
 m=pd.concat([ma,mb]);x=sp.vstack([a[:,common].X,b[:,common].X]).tocsr();assert m.donor.nunique()==55 and m.index.is_unique
 mmwrite(O/'counts.mtx',x.T);m.to_csv(O/'cell_metadata.tsv',sep='\t');pd.Series(common).to_csv(O/'genes.tsv',index=False,header=False,sep='\t')
 v={'status':'completed_input_with_partial_MTQC','cells':len(m),'genes':len(common),'donors':55,'cohort_cells':m.study.value_counts().to_dict(),'cohort_donors':m.groupby('study').donor.nunique().to_dict(),'F_donors':m.groupby('donor').F.first().value_counts().to_dict(),'excluded_ambiguous_G244832_features':int(dup.sum()),'QC_limitation':'G202379 MT absent and NA; G244832 MTQC existing. Shared gene intersection after each source QC; no zero-imputation of absent features.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
