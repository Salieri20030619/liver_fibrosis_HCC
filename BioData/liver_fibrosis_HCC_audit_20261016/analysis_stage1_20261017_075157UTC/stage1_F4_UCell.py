import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp,json,subprocess,traceback
from scipy.stats import rankdata,mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;S=P/'F4_official_sources_v1';O=P/'F4_UCell_Hallmark_v1';O.mkdir(exist_ok=False)
def bh(p):
 p=np.asarray(p);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
try:
 # Verify ported rank statistic against the actual author R function; no package substitution.
 helper=(S/'UCell_HelperFunctions.R').read_text();start=helper.index('u_stat <- function');end=helper.index('\n}',start)+2;function=helper[start:end]
 test=np.array([[1,2,4,9],[2.5,2.5,1,4],[4,3,2,1]],float);pd.DataFrame(test.T).to_csv(O/'test_ranks.tsv',sep='\t',index=False,header=False)
 rcode=function+'\nx<-read.table("'+str(O/'test_ranks.tsv')+'",header=FALSE); x<-cbind(rn=seq_len(nrow(x)),x); write.table(u_stat(x,c(1,3,-1),maxRank=4),file="'+str(O/'R_expected.tsv')+'",row.names=FALSE,col.names=FALSE)\n';(O/'author_function_test.R').write_text(rcode);r=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'author_function_test.R')],capture_output=True,text=True);assert r.returncode==0,r.stderr
 expect=np.loadtxt(O/'R_expected.tsv');actual=1-(np.minimum(test[:,[0,2]],4).sum(axis=1)+4-6)/(12-6);assert np.allclose(actual,expect);(O/'formula_validation.json').write_text(json.dumps({'matches_author_R_u_stat':True,'cases':3,'maxRank':4,'missing_gene_imputation':'maxRank','ties':'average ranks; no random tie-breaking','production_maxRank':1500},indent=2))
 sets={}
 for line in (S/'Hallmark_2025.1_Hs.gmt').read_text().splitlines():
  a=line.split('\t');sets[a[0]]=list(dict.fromkeys(a[2:]))
 assert len(sets)==50
 results=[];coverage=[]
 for g in ['GSE244832','GSE136103']:
  a=ad.read_h5ad(P/'F2_Roe_F3_F4_subclusters_v1'/(g+'_Myeloid')/'subcluster_checkpoint.h5ad');obs=a.obs.copy();eligible=obs.disease_comparison_eligible.astype(str).str.lower().eq('true');obs=obs[eligible].copy();obs['group']=np.where(obs.F_stage_from_Figure4C.astype(str).isin(['3','4','3.0','4.0']),'F3_F4','F0') if g=='GSE244832' else obs['disease status'];X=a.raw.X[eligible.to_numpy(),:].tocsr().copy();X.data=np.expm1(X.data);symbols=a.raw.var.gene_symbol_original.astype(str);unique=pd.Index(symbols.unique());codes=unique.get_indexer(symbols);combine=sp.csr_matrix((np.ones(len(codes)),(np.arange(len(codes)),codes)),shape=(len(codes),len(unique)));X=(X@combine).tocsr();maxrank=min(1500,len(unique));score=np.zeros((X.shape[0],50));indices=[]
  for name,genes in sets.items():
   ix=unique.get_indexer(genes);indices.append(ix);coverage.append({'cohort':g,'set':name,'genes':len(ix),'present':int((ix>=0).sum()),'missing':' ;'.join(np.array(genes)[ix<0])})
  for start in range(0,X.shape[0],128):
   ranks=np.minimum(rankdata(-X[start:start+128].toarray(),axis=1,method='average'),maxrank)
   for j,ix in enumerate(indices):
    n=len(ix);s=ranks[:,ix[ix>=0]].sum(axis=1)+int((ix<0).sum())*maxrank;score[start:start+128,j]=1-(s-n*(n+1)/2)/(n*maxrank-n*(n+1)/2)
  assert np.isfinite(score).all() and score.min()>-1e-9 and score.max()<1+1e-9
  scores=pd.DataFrame(score,index=obs.index,columns=list(sets));scores.to_csv(O/(g+'_cell_Hallmark_UCell.tsv.gz'),sep='\t');full=obs[['donor','sorting','group','subcluster_preliminary']].join(scores);donor=full.groupby(['donor','sorting','group'],observed=True)[list(sets)].mean().reset_index();donor.to_csv(O/(g+'_donor_Hallmark_UCell.tsv'),sep='\t',index=False);cluster=full.groupby(['donor','sorting','group','subcluster_preliminary'],observed=True)[list(sets)].mean().reset_index();cluster.to_csv(O/(g+'_donor_subcluster_Hallmark_UCell.tsv'),sep='\t',index=False)
  for sorting,ds in donor.groupby('sorting',observed=True):
   key=g+'_'+str(sorting).replace('+','positive').replace('-','negative');ref='F0' if g=='GSE244832' else 'healthy';case='F3_F4' if g=='GSE244832' else 'fibrotic';xx=ds[ds.group.eq(case)];yy=ds[ds.group.eq(ref)];rows=[]
   for name in sets:
    v=mannwhitneyu(xx[name],yy[name],alternative='two-sided',method='asymptotic');rows.append({'set':name,'case_donors':len(xx),'reference_donors':len(yy),'case_mean':xx[name].mean(),'reference_mean':yy[name].mean(),'mean_difference':xx[name].mean()-yy[name].mean(),'p_value':v.pvalue})
   t=pd.DataFrame(rows);t['BH_q']=bh(t.p_value);t.to_csv(O/(key+'_donor_pathway_tests.tsv'),sep='\t',index=False)
   cs=cluster[cluster.sorting.eq(sorting)];heat=cs.groupby(['group','subcluster_preliminary'],observed=True)[list(sets)].mean().T;z=heat.sub(heat.mean(axis=1),axis=0).div(heat.std(axis=1).replace(0,np.nan),axis=0).fillna(0);fig,ax=plt.subplots(figsize=(max(9,len(heat.columns)*.4),14));im=ax.imshow(z,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2);ax.set_yticks(range(50),[s.replace('HALLMARK_','') for s in sets],fontsize=7);ax.set_xticks(range(len(heat.columns)),[str(i[0])+':'+str(i[1]) for i in heat.columns],rotation=90,fontsize=7);ax.set_title(key+' | donor-balanced Hallmark UCell z-scores');fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(O/(key+'_F4DF_Hallmark_heatmap.png'),dpi=150);plt.close(fig)
   panel=['HALLMARK_GLYCOLYSIS','HALLMARK_HYPOXIA','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_OXIDATIVE_PHOSPHORYLATION'];fig,axes=plt.subplots(2,3,figsize=(13,7))
   for ax,name in zip(axes.ravel(),panel):
    ax.boxplot([yy[name],xx[name]],tick_labels=[ref,case]);ax.scatter(np.ones(len(yy)),yy[name],s=20);ax.scatter(np.full(len(xx),2),xx[name],s=20);tr=t[t.set.eq(name)].iloc[0];ax.set_title(name.replace('HALLMARK_','')+'\np='+str(round(tr.p_value,4))+' q='+str(round(tr.BH_q,4)),fontsize=9);ax.set_ylabel('Mean UCell per donor')
   fig.suptitle(key+' | individual points are donors');fig.tight_layout();fig.savefig(O/(key+'_F4E_donor_scores.png'),dpi=170);plt.close(fig);results.append({'comparison':key,'donors_case_reference':[len(xx),len(yy)],'pathways':50,'nominal_p_below_05':int((t.p_value<.05).sum()),'BH_q_below_05':int((t.BH_q<.05).sum()),'top5':t.sort_values('p_value').head(5).to_dict('records')})
  del a,X
 pd.DataFrame(coverage).to_csv(O/'gene_set_coverage.tsv',sep='\t',index=False);(O/'results.json').write_text(json.dumps(results,indent=2));(O/'method_scope.json').write_text(json.dumps({'method':'UCell rank-statistic Python port verified against downloaded author R u_stat; average ties, maxRank1500, missing genes imputed at maxRank. Exact formula, not claim execution of installed R UCell package.','version':'Official current helper source frozen by SHA256; Hallmark MSigDB2025.1.Hs. Article scoring version not specified.','aggregation':'sum duplicate symbols before rank; mean cell scores per donor/sorting, then equally weighted donor means within cluster for heatmap; no new cell filtering. Existing conflict exclusions retained.','panels':'F4D/F adapted Hallmark heatmaps and F4E donor score plots. Original custom M1/M2, lipid mediator, TCA sets not recovered; Hallmark OXPHOS is labelled OXPHOS, not substituted for TCA. F4G irGSEA remains unexecuted.','statistics':'two-sided donor rank-sum and BH over50 Hallmarks per cohort/sorting; no pooled-study test; no selected DTNA subgroup or mechanism inferred.'},indent=2));print(json.dumps({'status':'complete','comparisons':len(results)}),flush=True)
except Exception:(O/'error.txt').write_text(traceback.format_exc());raise
