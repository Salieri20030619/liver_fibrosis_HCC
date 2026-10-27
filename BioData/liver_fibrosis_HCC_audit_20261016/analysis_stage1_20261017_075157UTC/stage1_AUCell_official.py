import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,scipy.sparse as sp,scipy.io,json,subprocess,traceback
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'F4_AUCell_official_v1';O.mkdir(exist_ok=False);results=[]
for g in ['GSE244832','GSE136103']:
 try:
  d=O/g;d.mkdir();a=ad.read_h5ad(P/'F2_Roe_F3_F4_subclusters_v1'/(g+'_Myeloid')/'subcluster_checkpoint.h5ad');obs=a.obs.copy();keep=obs.disease_comparison_eligible.astype(str).str.lower().eq('true');obs=obs[keep].copy();obs['group']=np.where(obs.F_stage_from_Figure4C.astype(str).isin(['3','4','3.0','4.0']),'F3_F4','F0') if g=='GSE244832' else obs['disease status'];X=a.raw.X[keep.to_numpy(),:].tocsr().copy();X.data=np.expm1(X.data);symbols=a.raw.var.gene_symbol_original.astype(str);uni=pd.Index(symbols.unique());ix=uni.get_indexer(symbols);X=X@sp.csr_matrix((np.ones(len(ix)),(np.arange(len(ix)),ix)),shape=(len(ix),len(uni)));scipy.io.mmwrite(str(d/'expression.mtx'),X.T);(d/'genes.txt').write_text('\n'.join(uni)+'\n');(d/'cells.txt').write_text('\n'.join(obs.index)+'\n');obs.to_csv(d/'cell_metadata.tsv',sep='\t');del a,X
  r=f""".libPaths(c('{P/'trajectory_dependencies_v1/R_project_library_v1'}',.Library));library(AUCell);library(Matrix);set.seed(0)
d<-'{d}';x<-as(readMM(file.path(d,'expression.mtx')),'CsparseMatrix');rownames(x)<-readLines(file.path(d,'genes.txt'));colnames(x)<-readLines(file.path(d,'cells.txt'))
gmt<-strsplit(readLines('{P/'F4_official_sources_v1/Hallmark_2025.1_Hs.gmt'}'),'\\t');sets<-lapply(gmt,function(z)unique(z[-c(1,2)]));names(sets)<-sapply(gmt,function(z)z[1]);rankings<-AUCell_buildRankings(x,plotStats=FALSE,verbose=FALSE);auc<-AUCell_calcAUC(sets,rankings,aucMaxRank=ceiling(.05*nrow(x)),normAUC=TRUE,verbose=FALSE);write.table(t(getAUC(auc)),file.path(d,'official_AUCell_scores.tsv'),sep='\\t',quote=FALSE,col.names=NA);writeLines(capture.output(sessionInfo()),file.path(d,'sessionInfo.txt'))
""";(d/'run.R').write_text(r)
  with open(d/'R.log','w') as log:ret=subprocess.run(['/usr/bin/Rscript','--vanilla',str(d/'run.R')],stdout=log,stderr=subprocess.STDOUT)
  assert ret==0,('AUCell_R_failed',ret);scores=pd.read_csv(d/'official_AUCell_scores.tsv',sep='\t',index_col=0);assert scores.index.tolist()==obs.index.tolist() and np.isfinite(scores).all().all();full=obs[['donor','sorting','group']].join(scores);donor=full.groupby(['donor','sorting','group'],observed=True)[list(scores.columns)].mean().reset_index();donor.to_csv(d/'donor_AUCell_scores.tsv',sep='\t',index=False);u=pd.read_csv(P/'F4_UCell_Hallmark_v1'/(g+'_donor_Hallmark_UCell.tsv'),sep='\t');joined=donor.merge(u,on=['donor','sorting','group'],suffixes=('_AUC','_UCell'),validate='one_to_one');cor=[{'pathway':name,'Spearman_AUCell_vs_UCell':joined[name+'_AUC'].corr(joined[name+'_UCell'],method='spearman')} for name in scores.columns];pd.DataFrame(cor).to_csv(d/'AUCell_UCell_donor_concordance.tsv',sep='\t',index=False)
  for sort,z in donor.groupby('sorting',observed=True):
   groups=sorted(z.group.unique());rows=[]
   for name in scores.columns:
    xx=z[z.group.eq(groups[0])][name];yy=z[z.group.eq(groups[1])][name];rows.append({'pathway':name,'group1':groups[0],'group2':groups[1],'n1':len(xx),'n2':len(yy),'mean_difference':xx.mean()-yy.mean(),'p_value':mannwhitneyu(xx,yy,alternative='two-sided',method='asymptotic').pvalue})
   t=pd.DataFrame(rows);p=t.p_value.to_numpy();order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);t['BH_q']=q;key=str(sort).replace('+','positive').replace('-','negative');t.to_csv(d/(key+'_pathway_tests.tsv'),sep='\t',index=False);show=['HALLMARK_GLYCOLYSIS','HALLMARK_HYPOXIA','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_COAGULATION'];fig,axs=plt.subplots(1,4,figsize=(14,4))
   for ax,name in zip(axs,show):
    vals=[z[z.group.eq(h)][name] for h in groups];ax.boxplot(vals,tick_labels=groups)
    for i,v in enumerate(vals):ax.scatter(np.full(len(v),i+1),v,s=16)
    ax.set_title(name.replace('HALLMARK_',''),fontsize=8);ax.set_ylabel('Mean official AUCell per donor')
   fig.suptitle(g+' '+str(sort)+' | AUCell companion, not irGSEA integration');fig.tight_layout();fig.savefig(d/(key+'_F4E_AUCell_companion.png'),dpi=160);plt.close(fig);results.append({'cohort':g,'sorting':str(sort),'BH_q_below05':int((q<.05).sum()),'donors':len(z)})
  (O/'state.json').write_text(json.dumps(results));print(json.dumps(results),flush=True)
 except Exception:(d/'error.txt').write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
(O/'scope.json').write_text(json.dumps({'method':'Official AUCell1.32.0; seed0; aucMaxRank ceiling5% genes; normAUC TRUE. Default rank cutoff is algorithm parameter, no new cell filter.','panels':'F4E companion; full irGSEA F4G remains separate','input':'Same eligible myeloid cells, duplicate symbols summed before ranks; public Hallmark2025.1.Hs, not original unavailable custom sets','statistical_unit':'donor x sorting; BH across50 pathways per comparison; disease labels and D9/D10 exclusions retained'},indent=2))
