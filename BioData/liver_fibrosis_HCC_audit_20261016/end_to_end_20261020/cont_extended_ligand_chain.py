from pathlib import Path
import time,json,shutil,numpy as np,pandas as pd
from scipy.stats import spearmanr,mannwhitneyu
from statsmodels.stats.multitest import multipletests
import gseapy as gp
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';CC=P/'continuation_CellChat_expanded_v1';TF=P/'continuation_SCENIC_candidates_v1';O=P/'continuation_ligand_to_GSEA_v1';O.mkdir(exist_ok=False)
for i in range(1080):
 if (CC/'validation.json').exists() and (TF/'validation.json').exists():break
 if (P/'continuation_SCENIC_expanded_v1/error.txt').exists():raise RuntimeError('SCENIC upstream failed; downstreamnotrun')
 time.sleep(10)
else:raise RuntimeError('Upstream pending3hours')
cc=pd.read_csv(CC/'all_record_state_specific_communications.tsv.gz',sep='\t');cc=cc[(cc.F.astype(float)>=3)&cc.source.str.startswith('Myeloid_')&cc.target.eq('Stromal')&cc.pval.le(.05)];rank=cc.groupby('ligand').agg(source_records=('source_record','nunique'),mean_probability=('prob','mean')).sort_values(['source_records','mean_probability'],ascending=False);rank.to_csv(O/'advanced_myeloid_to_stromal_ligand_ranking.tsv',sep='\t')
if not len(rank):
 (O/'validation.json').write_text(json.dumps({'status':'completed_negative_no_eligible_ligand','scope':'No significantadvancedF3/F4myeloid-state-to-stromal interactions; no forcedligand or downstreamGSEA'},indent=2));raise SystemExit(0)
x=pd.read_csv(D/'GSE162694_raw_counts.csv.gz',index_col=0);m=pd.read_csv(D/'expression_metadata_alignment.tsv',sep='\t',keep_default_na=False).set_index('expression_id').loc[x.columns];m=m[m['fibrosis stage'].isin(['0','1','2','3','4'])];x=x[m.index];h=pd.read_csv(D/'HGNC_one_to_one_mapping_v1.tsv',sep='\t').set_index('ensembl_gene_id').symbol;norm=np.log2(x.div(x.sum(0),axis=1)*1e6+1);names=pd.Series(x.index,index=x.index).map(h);norm=norm.loc[names.notna()];norm.index=names.dropna();assert norm.index.is_unique
available=[g for g in rank.index if g in norm.index];assert available,'No selectedligand hasunique publicbulkfeature';ligand=available[0];threshold=float(norm.loc[ligand].quantile(.75));hi=norm.loc[ligand].ge(threshold);assert min(hi.sum(),(~hi).sum())>=10;design=m.copy();design['ligand_expression']=norm.loc[ligand];design['ligand_high']=hi;design.to_csv(O/'ligand_sample_design.tsv',sep='\t');corr=[]
for gene in pd.read_csv(TF/'F5G_candidates.tsv',sep='\t').TF.unique():
 if gene in norm.index:
  rho,p=spearmanr(norm.loc[gene],norm.loc[ligand]);corr.append({'TF':gene,'ligand':ligand,'rho':rho,'p':p})
c=pd.DataFrame(corr,columns=['TF','ligand','rho','p']);c['BH_q']=multipletests(c.p,method='fdr_bh')[1] if len(c) else [];c.to_csv(O/'new_network_TF_ligand_correlations.tsv',sep='\t',index=False)
old=P/'F6_ligand_to_F7_GSEA';reuse=json.loads((old/'validation.json').read_text())['ligand']==ligand
if reuse:
 for f in ['ligand_high_vs_low_DE.tsv','full_ranked_genes.tsv','Hallmark_GSEA.tsv']:shutil.copyfile(old/f,O/f)
else:
 aa=norm.loc[:,hi].to_numpy();bb=norm.loc[:,~hi].to_numpy();test=mannwhitneyu(aa,bb,axis=1,alternative='two-sided',method='asymptotic');delta=aa.mean(1)-bb.mean(1);de=pd.DataFrame({'gene':norm.index,'mean_log2CPM_difference':delta,'p':test.pvalue});de['BH_q']=multipletests(de.p,method='fdr_bh')[1];de['rank_score']=np.sign(delta)*-np.log10(np.maximum(de.p,1e-300));de.to_csv(O/'ligand_high_vs_low_DE.tsv',sep='\t',index=False);rnk=pd.read_csv(O/'ligand_high_vs_low_DE.tsv',sep='\t')[['gene','rank_score']].sort_values(['rank_score','gene'],ascending=[False,True]);rnk.to_csv(O/'full_ranked_genes.tsv',sep='\t',index=False);res=gp.prerank(rnk=rnk,gene_sets=str(P/'Hallmark_2025.1_Hs.gmt'),min_size=15,max_size=500,permutation_num=1000,threads=4,seed=20261020,outdir=None,no_plot=True);res.res2d.to_csv(O/'Hallmark_GSEA.tsv',sep='\t',index=False)
v={'status':'completed','ligand':ligand,'advanced_source_records_contributing':cc.source_record.nunique(),'states_contributing':cc.source.nunique(),'ligand_high':int(hi.sum()),'ligand_low':int((~hi).sum()),'new_network_TFs_tested':len(c),'DE_GSEA_reused_identical_ligand_and_same_bulk_design':reuse,'selection':'Highest source-recordrecurrent significantadvancedF3/F4 myeloidstate-to-stromal ligand withuniqueG162694feature; noPLGpreset. NewmotifnetworkTFcandidatefile physicallyconsumed.','limitations':'AmbienthepaticRNA remains unresolvedwithoutemptydroplets; secretedhepaticligandsinmyeloid areassociations requiringvalidation. StromalnotpureaHSC. Crosssectionalbulkcomposition/Fconfounding; notincidentHCC or causality.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
