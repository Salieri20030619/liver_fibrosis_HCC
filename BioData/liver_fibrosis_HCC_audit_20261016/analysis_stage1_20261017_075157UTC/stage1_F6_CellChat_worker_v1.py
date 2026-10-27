import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,subprocess,hashlib,datetime,traceback
import numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
from scipy.io import mmwrite
P=Path(__file__).parent;O=P/'F6_CellChat_GSE244832_v1';O.mkdir(exist_ok=False)
L=O/'R_project_library';L.mkdir();N=P/'dependency_irgsea_repair_v3/native_dependencies';env=os.environ.copy();env['LD_PRELOAD']=str(N/'lib/libstdc++.so.6')+':'+str(N/'lib/libdeflate.so.0');env['LD_LIBRARY_PATH']=str(N/'lib');env['PKG_CONFIG_PATH']=str(N/'lib/pkgconfig');env['MAKEFLAGS']='-j4';env['USE_BUNDLED_LIBUV']='1'
libs=f".libPaths(c('{L}','{P/'irGSEA_dependencies_v1/R_project_library_v1'}','{P/'trajectory_dependencies_v1/R_project_library_v1'}',.Library));"
def state(stage,**kw):
 v={'stage':stage,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw};(O/'state.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
try:
 state('preparing donor inputs')
 m=pd.read_csv(P/'coarse_celltype_review_v2/GSE244832_reviewed_cell_metadata.tsv.gz',sep='\t',index_col=0,keep_default_na=False)
 a=ad.read_h5ad(P/'GSE244832/qc_pass_counts.h5ad');assert m.index.tolist()==a.obs_names.tolist()
 eligible=m.disease_comparison_eligible.astype(str).str.lower().eq('true');assert m.loc[eligible,'donor'].nunique()==16
 keep=eligible & m.coarse_review_status.str.startswith('supported_coarse');m=m.loc[keep].copy();X=a.X[keep.to_numpy()].tocsr().astype(np.float64);symbols=a.var_names.astype(str);assert symbols.is_unique
 assert np.isfinite(X.data).all() and (X.data>=0).all() and np.allclose(X.data,np.rint(X.data))
 pd.DataFrame({'gene':symbols}).to_csv(O/'gene_ids.tsv',sep='\t',index=False)
 targets=['Myeloid_macrophage_like','Mesenchymal_stromal','Mesenchymal_mural'];m['communication_label']=m.coarse_label_reviewed
 m[['donor','library','condition_Excel_original','F_stage_from_Figure4C','coarse_label_reviewed']].to_csv(O/'eligible_cell_metadata.tsv.gz',sep='\t')
 rules={'cohort':'GSE244832','clinical_exclusions':['D9','D10'],'identity_resolution':targets,'identity_source':'coarse_celltype_review_v2; canonical multi-marker and donor review, no forced HSC/DTNA subtype','CellChat':'Official fixed source commit75253cd0c9e68410e6e721a6d3a0419a1d7e358f; human full bundled DB; log1pCP10k from full gene counts per cell; each donor separately; default probability parameters, min.cells10','comparison':'Descriptive donor-specific networks; cell-level permutation p is not a disease-level p; NORMAL/MASL/MASH retained separately; no pairing across people','signature':'Before testing bulk outcomes: target vs all other supported broad cells within donor, minimum20 target and20 rest cells; mean log2(CPM+1) pseudobulk difference ranking; positive in>=75% eligible donors, target mean detection>=25%, at least8 eligible donors; top50 ranked by mean effect (gene name tiebreak); preserve full rankings. No disease significance or original axis used.','limitations':['Ambient hepatic RNA/heterogeneity remains unresolved','Stromal/mural are broad compartments, not pure aHSC states','No probability proves physical signaling or causal mechanism','Original CellChat version/database and exact signature method not specified; reconstructed choices disclosed']}
 (O/'analysis_rules_frozen.json').write_text(json.dumps(rules,indent=2));counts=m.groupby(['donor','communication_label']).size().unstack(fill_value=0);counts.to_csv(O/'donor_celltype_counts.tsv',sep='\t')
 S=O/'signature_provenance';S.mkdir();siglines=[];sigstats=[]
 for target in targets:
  effects=[];detections=[];donors=[]
  for donor,z in m.groupby('donor'):
   pos=m.index.get_indexer(z.index);tt=z.coarse_label_reviewed.eq(target).to_numpy();rr=~tt
   if tt.sum()<20 or rr.sum()<20:continue
   tx=X[pos[tt]];rx=X[pos[rr]];u=np.asarray(tx.sum(axis=0)).ravel();v=np.asarray(rx.sum(axis=0)).ravel();effects.append(np.log2(u/u.sum()*1e6+1)-np.log2(v/v.sum()*1e6+1));detections.append(np.asarray((tx>0).mean(axis=0)).ravel());donors.append(donor)
  d=np.asarray(effects);det=np.asarray(detections);tab=pd.DataFrame({'gene':symbols,'mean_log2CPM_difference':d.mean(axis=0),'positive_donor_fraction':(d>0).mean(axis=0),'mean_target_detection':det.mean(axis=0),'eligible_donors':len(donors)})
  for i,donor in enumerate(donors):tab['effect_'+donor]=d[i]
  tab['passes_frozen_rule']=(tab.positive_donor_fraction>=.75)&(tab.mean_target_detection>=.25)&(len(donors)>=8)
  tab=tab.sort_values(['mean_log2CPM_difference','gene'],ascending=[False,True]);selected=tab[tab.passes_frozen_rule].head(50).gene.tolist();tab['selected']=tab.gene.isin(selected);tab.to_csv(S/(target+'_all_gene_ranking.tsv.gz'),sep='\t',index=False);siglines.append(target+'\tGSE244832_donor_consistent_broad_identity\t'+'\t'.join(selected));sigstats.append({'signature':target,'donors':donors,'selected_genes':selected,'eligible_genes':int(tab.passes_frozen_rule.sum())})
 (S/'broad_signatures.gmt').write_text('\n'.join(siglines)+'\n');(S/'signature_audit.json').write_text(json.dumps(sigstats,indent=2))
 inputs=[]
 for donor,z in m.groupby('donor'):
  labcounts=z.coarse_label_reviewed.value_counts();labels=[t for t in targets if labcounts.get(t,0)>=10];sub=z[z.coarse_label_reviewed.isin(labels)]
  if len(labels)<2:inputs.append({'donor':donor,'status':'skipped_less_than_two_eligible_labels','labels':labels});continue
  ix=m.index.get_indexer(sub.index);xx=X[ix].copy();tot=np.asarray(xx.sum(axis=1)).ravel();xx=sp.diags(1e4/tot)@xx;xx.data=np.log1p(xx.data);q=O/donor;q.mkdir();mmwrite(str(q/'normalized_expression.mtx'),xx.T);pd.DataFrame({'cell':sub.index,'labels':sub.coarse_label_reviewed.to_numpy()}).to_csv(q/'cell_metadata.tsv',sep='\t',index=False);inputs.append({'donor':donor,'status':'ready','cells':len(sub),'labels':labels,'condition':sub.condition_Excel_original.iloc[0],'F_stage':int(sub.F_stage_from_Figure4C.iloc[0])})
 (O/'input_manifest.json').write_text(json.dumps(inputs,indent=2));state('inputs and signatures ready',donors=len(inputs))
 tar=P/'CellChat_75253cd0c9e68410e6e721a6d3a0419a1d7e358f.tar.gz';assert hashlib.sha256(tar.read_bytes()).hexdigest()=='52a90ce9f121acd58c9e5b2542956de0f3d81ac9b10906454847946f289e0739'
 r=libs+f"options(repos=c(CRAN='https://cloud.r-project.org'),timeout=300);pk<-c('NMF','ggalluvial','svglite','sna','ggnetwork','collapse');missing<-pk[!sapply(pk,requireNamespace,quietly=TRUE)];if(length(missing))install.packages(missing,lib='{L}',Ncpus=4,dependencies=NA);install.packages('{tar}',repos=NULL,type='source',lib='{L}');stopifnot(requireNamespace('CellChat',quietly=TRUE));writeLines(capture.output(sessionInfo()),'{O/'installation_sessionInfo.txt'}')"
 (O/'install.R').write_text(r);state('installing CellChat in new isolated project library')
 with open(O/'install.log','w') as log:ret=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'install.R')],stdout=log,stderr=subprocess.STDOUT,env=env)
 if ret.returncode:raise RuntimeError('CellChat installation failed; see install.log')
 state('running donor-specific CellChat')
 r=libs+f"""library(CellChat);library(Matrix);library(future);plan('sequential');o<-'{O}';genes<-read.delim(file.path(o,'gene_ids.tsv'))$gene;data(CellChatDB.human);saveRDS(CellChatDB.human,file.path(o,'CellChatDB.human.rds'));write.table(CellChatDB.human$interaction,file.path(o,'database_interactions.tsv'),sep='\t',quote=FALSE,row.names=FALSE);writeLines(capture.output(sessionInfo()),file.path(o,'sessionInfo.txt'))
 dirs<-list.dirs(o,recursive=FALSE,full.names=TRUE);dirs<-dirs[file.exists(file.path(dirs,'normalized_expression.mtx'))]
 for(d in dirs){{
  tryCatch({{
   x<-readMM(file.path(d,'normalized_expression.mtx'));meta<-read.delim(file.path(d,'cell_metadata.tsv'),row.names=1);rownames(x)<-genes;colnames(x)<-rownames(meta);stopifnot(identical(colnames(x),rownames(meta)));obj<-createCellChat(x,meta=meta,group.by='labels');obj@DB<-CellChatDB.human;obj<-subsetData(obj);obj<-identifyOverExpressedGenes(obj);obj<-identifyOverExpressedInteractions(obj);obj<-computeCommunProb(obj);obj<-filterCommunication(obj,min.cells=10);obj<-computeCommunProbPathway(obj);obj<-aggregateNet(obj);saveRDS(obj,file.path(d,'CellChat.rds'));tab<-subsetCommunication(obj);write.table(tab,file.path(d,'communications.tsv'),sep='\t',quote=FALSE,row.names=FALSE);write.table(obj@net$count,file.path(d,'network_count.tsv'),sep='\t',quote=FALSE,col.names=NA);write.table(obj@net$weight,file.path(d,'network_weight.tsv'),sep='\t',quote=FALSE,col.names=NA);writeLines('completed',file.path(d,'status.txt'))
  }},error=function(e){{writeLines(conditionMessage(e),file.path(d,'error.txt'))}})
 }}
 """;(O/'run_CellChat.R').write_text(r)
 with open(O/'CellChat.log','w') as log:ret=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run_CellChat.R')],stdout=log,stderr=subprocess.STDOUT,env=env)
 statuses={z['donor']:('completed' if (O/z['donor']/'status.txt').exists() else ((O/z['donor']/'error.txt').read_text() if (O/z['donor']/'error.txt').exists() else z['status'])) for z in inputs};state('finished',R_returncode=ret.returncode,donors=statuses)
except Exception:
 (O/'worker_error.txt').write_text(traceback.format_exc());state('failed',error=traceback.format_exc());raise
