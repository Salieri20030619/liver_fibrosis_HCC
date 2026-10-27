from pathlib import Path
import pandas as pd,numpy as np,json,gseapy as gp
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent; L=P/'continuation_ligand_to_GSEA_v1';O=P/'continuation_GSEA_rank_sensitivity_v1';O.mkdir(exist_ok=False)
d=pd.read_csv(L/'ligand_high_vs_low_DE.tsv',sep='\t');x=pd.read_csv(P/'F6_S4_S6_F8_gene_validation/GSE162694_log2CPM_full.tsv.gz',sep='\t',index_col=0).loc[d.gene];keep=(x>1).sum(axis=1).ge(np.ceil(x.shape[1]*.1));rows=[];outs=[]
for universe,dd in [('all',d),('expressed_CPM1_10percent',d.loc[keep.to_numpy()])]:
 dd=dd.copy();dd['BH_q_within_universe']=multipletests(dd.p,method='fdr_bh')[1];dd.to_csv(O/(universe+'_genes.tsv'),sep='\t',index=False)
 vals=np.sort(dd.rank_score.unique());gap=np.diff(vals).min();assert gap>0
 for seed in [20261020,20261021]:
  # Strictly bounded perturbation breaks only exact ties, never reverses distinct original scores.
  r=dd[['gene','rank_score']].copy();rng=np.random.default_rng(seed);r['rank_score']+=rng.uniform(-gap/4,gap/4,len(r));r=r.sort_values('rank_score',ascending=False);r.to_csv(O/(universe+'_'+str(seed)+'_ranks.tsv'),sep='\t',index=False)
  res=gp.prerank(rnk=r,gene_sets=str(P/'Hallmark_2025.1_Hs.gmt'),min_size=15,max_size=500,permutation_num=1000,threads=4,seed=20261020,outdir=None,no_plot=True).res2d
  res.to_csv(O/(universe+'_'+str(seed)+'_GSEA.tsv'),sep='\t',index=False);res['universe']=universe;res['tie_seed']=seed;outs.append(res)
  rows.append({'universe':universe,'seed':seed,'genes':len(dd),'duplicate_fraction_original':float(dd.rank_score.duplicated().mean()),'minimum_distinct_score_gap':gap,'BH_significant_DE':int((dd.BH_q_within_universe<.05).sum())})
allr=pd.concat(outs);allr.to_csv(O/'all_sensitivity_GSEA.tsv',sep='\t',index=False)
stab=allr.groupby('Term').agg(NES_min=('NES','min'),NES_max=('NES','max'),FDR_min=('FDR q-val','min'),FDR_max=('FDR q-val','max'));stab['same_NES_sign']=np.sign(stab.NES_min)==np.sign(stab.NES_max);stab['FDR_lt05_all_four']=stab.FDR_max.lt(.05);stab.to_csv(O/'pathway_stability.tsv',sep='\t')
v={'status':'completed','runs':rows,'pathways':len(stab),'same_direction_all_four':int(stab.same_NES_sign.sum()),'FDR_lt05_all_four':int(stab.FDR_lt05_all_four.sum()),'reason':'Mann Whitney rank statistics are discrete; 78.47 percent ties do not arise from all-zero genes. Original results retained. Two reproducible tie orders with fixed GSEA permutation seed and expression-filter sensitivity; no DEG-only ranking.','limits':'Sensitivity is exploratory and post-hoc motivated by warning; no independent validation, no causal or incident HCC inference.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
