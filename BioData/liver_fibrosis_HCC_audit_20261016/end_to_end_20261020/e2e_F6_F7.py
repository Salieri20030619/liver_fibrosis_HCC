from pathlib import Path
import json,traceback
import pandas as pd,numpy as np
from scipy.stats import spearmanr,mannwhitneyu
from statsmodels.stats.multitest import multipletests
import gseapy as gp
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'F6_ligand_to_F7_GSEA';O.mkdir(exist_ok=False)
try:
 cc=pd.read_csv(P/'all_donor_LR_results.tsv.gz',sep='\t');cc=cc[cc.donor.isin(['D'+str(i) for i in range(11,19)])&cc.source.eq('Myeloid_macrophage_like')&cc.target.eq('Mesenchymal_stromal')&cc.pval.le(.05)]
 rank=cc.groupby('ligand').agg(donors=('donor','nunique'),mean_probability=('prob','mean')).sort_values(['donors','mean_probability'],ascending=False);assert len(rank)>0;rank.to_csv(O/'source_ligand_ranking.tsv',sep='\t')
 x=pd.read_csv(D/'GSE162694_raw_counts.csv.gz',index_col=0);m=pd.read_csv(D/'expression_metadata_alignment.tsv',sep='\t',keep_default_na=False).set_index('expression_id').loc[x.columns];m=m[m['fibrosis stage'].isin(['0','1','2','3','4'])];x=x[m.index]
 h=pd.read_csv(D/'HGNC_one_to_one_mapping_v1.tsv',sep='\t').set_index('ensembl_gene_id').symbol;norm=np.log2(x.div(x.sum(axis=0),axis=1)*1e6+1);names=pd.Series(x.index,index=x.index).map(h);norm=norm.loc[names.notna()];norm.index=names.dropna();assert norm.index.is_unique
 ligand=next(v for v in rank.index if v in norm.index);threshold=float(norm.loc[ligand].quantile(.75));hi=norm.loc[ligand].ge(threshold);assert hi.sum()>=10 and (~hi).sum()>=10
 design=m.copy();design['ligand_expression']=norm.loc[ligand];design['ligand_high']=hi;design.to_csv(O/'ligand_sample_design.tsv',sep='\t')
 tf=pd.read_csv(P/'F5_to_F6_original_direction_v2/F5G_candidates.tsv',sep='\t').TF.unique();corr=[]
 for t in tf:
  if t in norm.index:
   rho,p=spearmanr(norm.loc[t],norm.loc[ligand]);corr.append({'TF':t,'ligand':ligand,'rho':rho,'p':p})
 c=pd.DataFrame(corr);c['BH_q']=multipletests(c.p,method='fdr_bh')[1];c.to_csv(O/'TF_ligand_correlations.tsv',sep='\t',index=False)
 a=norm.loc[:,hi].to_numpy();b=norm.loc[:,~hi].to_numpy();test=mannwhitneyu(a,b,axis=1,alternative='two-sided',method='asymptotic');delta=a.mean(axis=1)-b.mean(axis=1)
 de=pd.DataFrame({'gene':norm.index,'mean_log2CPM_difference':delta,'p':test.pvalue});de['BH_q']=multipletests(de.p,method='fdr_bh')[1];de['rank_score']=np.sign(delta)*-np.log10(np.maximum(de.p,1e-300));de.to_csv(O/'ligand_high_vs_low_DE.tsv',sep='\t',index=False)
 # GSEA consumes the saved complete ranked output, not only significant genes.
 saved=pd.read_csv(O/'ligand_high_vs_low_DE.tsv',sep='\t');rnk=saved[['gene','rank_score']].sort_values(['rank_score','gene'],ascending=[False,True]);rnk.to_csv(O/'full_ranked_genes.tsv',sep='\t',index=False)
 result=gp.prerank(rnk=rnk,gene_sets=str(P/'Hallmark_2025.1_Hs.gmt'),min_size=15,max_size=500,permutation_num=1000,threads=4,seed=20261020,outdir=None,no_plot=True,verbose=True);result.res2d.to_csv(O/'Hallmark_GSEA.tsv',sep='\t',index=False)
 v={'status':'completed','source_discovery_donors':8,'ligand':ligand,'ligand_high':int(hi.sum()),'ligand_low':int((~hi).sum()),'genes_ranked':len(rnk),'TF_correlations':len(c),'pathways_tested':len(result.res2d),'permutations':1000,'adaptations':'Highest recurrent discovered myeloid-to-stromal ligand, no PLG prespecified; GSE162694112F-labelled samples; 75th-percentile split as original; full signed-log-p Wilcoxon ranking reconstructed; Hallmark2025.1','limitations':'Broad stromal label not pure aHSC; cross-sectional expression associations can reflect composition/F; no causality or incident-HCC prediction. Spatial branch separate and not replaced.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
