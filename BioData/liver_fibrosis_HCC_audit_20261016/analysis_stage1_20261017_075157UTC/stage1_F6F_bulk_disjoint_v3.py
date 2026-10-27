from pathlib import Path
import os,json,subprocess,hashlib,itertools
import numpy as np,pandas as pd
from scipy.stats import spearmanr,false_discovery_control
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'F6F_GSE164760_signature_disjoint_v3';O.mkdir(exist_ok=False);source=P/'F6_CellChat_GSE244832_v1/signature_provenance/broad_signatures.gmt';assert source.exists()
gmt=source.read_text();(O/'source_signatures.gmt').write_text(gmt);sets={r[0]:r[2:] for line in gmt.splitlines() if (r:=line.split('\t'))};assert len(sets)==3 and all(len(v)>=10 for v in sets.values())
original_sets={k:list(v) for k,v in sets.items()}
shared={g for v in sets.values() for g in v if sum(g in other for other in sets.values())>1}
sets={k:[g for g in v if g not in shared] for k,v in sets.items()}
(O/'removed_shared_genes.json').write_text(json.dumps({'removed':sorted(shared),'original_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'purpose':'Predefined overlap sensitivity; remove shared genes from every signature, preserve other parameters and all tissue groups'},indent=2))
(O/'source_signatures.gmt').write_text('\n'.join(k+'\tGSE244832_disjoint_sensitivity\t'+'\t'.join(v) for k,v in sets.items())+'\n')
overlap=[]
for a,b in itertools.combinations(sets,2):overlap.append({'signature1':a,'signature2':b,'shared_genes':sorted(set(sets[a])&set(sets[b]))})
(O/'signature_overlap.json').write_text(json.dumps(overlap,indent=2))
meta=pd.read_csv(P/'GSE164760_design_eligibility_v1/sample_design.tsv',sep='\t',keep_default_na=False).set_index('sample_id');x=pd.read_csv(P/'GSE164760_gene_expression_median_v1.tsv.gz',sep='\t',index_col=0);assert x.columns.tolist()==meta.index.tolist()
coverage=[]
for k,v in sets.items():coverage.append({'signature':k,'source_genes':len(v),'matched_genes':len(set(v)&set(x.index)),'missing':sorted(set(v)-set(x.index))})
assert all(v['matched_genes']>=10 for v in coverage);(O/'signature_coverage.json').write_text(json.dumps(coverage,indent=2))
rules={'panel':'F6F exploratory adaptation, not exact original abundance estimation','source':'GSE244832 donor-consistent broad signatures frozen without bulk outcomes or disease p-value selection','statistical_plan':'All3 signature pairs within each of5 tissue groups; Spearman rho as descriptive association only. NO p/q inferred because individual patient independence and tumor/adjacent linkage unresolved. No pooled across-group correlation.','limitations':['Scores are expression signatures, not deconvolved cell fractions','No same-patient assumption for tumor/adjacent','No longitudinal fibrosis-to-HCC inference','Mixed/ambient hepatic signal may affect source signatures','No candidate axis assumed','Single HCC bulk cohort, not nine independently replicated cohorts'],'source_signature_sha256':hashlib.sha256(source.read_bytes()).hexdigest()};(O/'analysis_rules.json').write_text(json.dumps(rules,indent=2))
r=f""".libPaths(c('{P/'irGSEA_dependencies_v1/R_project_library_v1'}','{P/'trajectory_dependencies_v1/R_project_library_v1'}',.Library));library(GSVA);library(BiocParallel);x<-as.matrix(read.delim(gzfile('{P/'GSE164760_gene_expression_median_v1.tsv.gz'}'),row.names=1,check.names=FALSE));x<-log2(x);g<-strsplit(readLines('{O/'source_signatures.gmt'}'),'\t');s<-lapply(g,function(z)unique(z[-c(1,2)]));names(s)<-sapply(g,function(z)z[1]);z<-gsva(ssgseaParam(x,s,minSize=1,maxSize=Inf,normalize=TRUE),BPPARAM=SerialParam());write.table(z,'{O/'official_ssGSEA_scores.tsv'}',sep='\t',quote=FALSE,col.names=NA);writeLines(capture.output(sessionInfo()),'{O/'sessionInfo.txt'}')""";(O/'run.R').write_text(r);n=P/'dependency_irgsea_repair_v3/native_dependencies';env=os.environ.copy();env['LD_PRELOAD']=str(n/'lib/libstdc++.so.6')+':'+str(n/'lib/libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n/'lib');env['OMP_NUM_THREADS']='4';env['OPENBLAS_NUM_THREADS']='4'
with open(O/'R.log','w') as log:ret=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run.R')],stdout=log,stderr=subprocess.STDOUT,env=env)
assert ret.returncode==0;sc=pd.read_csv(O/'official_ssGSEA_scores.tsv',sep='\t',index_col=0);assert sc.columns.tolist()==meta.index.tolist() and sc.shape==(3,170) and np.isfinite(sc).all().all();sc.T.join(meta).to_csv(O/'scores_with_metadata.tsv',sep='\t');rows=[]
fig,axes=plt.subplots(1,5,figsize=(15,4));groups=sorted(meta.group.unique())
for ax,group in zip(axes,groups):
 ids=meta.index[meta.group.eq(group)];z=sc[ids];mat=np.eye(3)
 for a,b in itertools.combinations(range(3),2):
  rho=float(spearmanr(z.iloc[a],z.iloc[b]).statistic);assert np.isclose(rho,np.corrcoef(z.iloc[a].rank(),z.iloc[b].rank())[0,1]);mat[a,b]=mat[b,a]=rho;rows.append({'group':group,'samples':len(ids),'signature1':sc.index[a],'signature2':sc.index[b],'Spearman_rho':rho,'inference':'descriptive only; no patient-independent p/q'})
 im=ax.imshow(mat,cmap='RdBu_r',vmin=-1,vmax=1);ax.set_xticks(range(3),['Myeloid','Stromal','Mural'],rotation=60);ax.set_yticks(range(3),['Myeloid','Stromal','Mural']);ax.set_title(group.replace('Non-tumoral NASH liver adjacent to HCC','HCC adjacent').replace(' liver','')+'\nn='+str(len(ids)),fontsize=9)
 for i in range(3):
  for j in range(3):ax.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',fontsize=8)
fig.suptitle('GSE164760: within-group broad signature correlations; descriptive only');fig.tight_layout(rect=[0,0,1,.92]);fig.savefig(O/'within_group_signature_correlations.png',dpi=180);plt.close(fig);pd.DataFrame(rows).to_csv(O/'within_group_correlations.tsv',sep='\t',index=False)
v={'status':'completed','samples':170,'scores':3,'groups':5,'correlations':15,'all_IDs_exact':True,'finite_scores':True,'independent_rank_correlation_validation':True,'p_values_computed':False,'signature_overlap':overlap};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
