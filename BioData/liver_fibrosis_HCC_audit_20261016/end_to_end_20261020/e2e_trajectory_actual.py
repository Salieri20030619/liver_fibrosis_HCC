from pathlib import Path
import subprocess,os,json
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'F3_F5_ordered_trajectory';O.mkdir(exist_ok=False)
r=r'''.libPaths(c('LIB','OLD/trajectory_dependencies_v1/R_project_library_v1','OLD/irGSEA_dependencies_v1/R_project_library_v1',.libPaths()));library(monocle);library(CytoTRACE);library(Matrix)
for(branch in c('F5','F3')){od<-paste0('OUTPUT/',branch);dir.create(od);tryCatch({
 src<-if(branch=='F5')'OLD/F5_myeloid_potential_topology_v1' else 'OLD/F3_DDRTree_unrooted_v1/GSE244832'
 cds<-readRDS(file.path(src,'DDRTree_unrooted.rds'))
 if(branch=='F3'){
  keep<-tolower(as.character(pData(cds)$disease_comparison_eligible))=='true';cds<-cds[,keep];set.seed(0);cds<-reduceDimension(cds,max_components=2,reduction_method='DDRTree',norm_method='log',pseudo_expr=1,relative_expr=TRUE,verbose=FALSE)
  z<-CytoTRACE::CytoTRACE(as.matrix(exprs(cds)),enableFast=TRUE,ncores=2,subsamplesize=1000);cyto<-z$CytoTRACE;saveRDS(z,file.path(od,'CytoTRACE_core.rds'))
 }else{z<-read.delim(file.path(src,'cell_scores.tsv'),row.names=1);cyto<-setNames(z$CytoTRACE,rownames(z))}
 cds<-orderCells(cds);pData(cds)$CytoTRACE<-cyto[rownames(pData(cds))];stopifnot(all(is.finite(pData(cds)$CytoTRACE)))
 dm<-aggregate(CytoTRACE~State+donor,data=pData(cds),mean);state<-aggregate(CytoTRACE~State,data=dm,mean);root<-state$State[which.max(state$CytoTRACE)];write.table(state,file.path(od,'root_selection_equal_donor.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
 cds<-orderCells(cds,root_state=root);write.table(pData(cds),file.path(od,'ordered_cell_metadata.tsv'),sep='\t',quote=FALSE,col.names=NA);write.table(table(pData(cds)$State,pData(cds)$donor),file.path(od,'state_donor_counts.tsv'),sep='\t',quote=FALSE)
 cds<-estimateDispersions(cds);genes<-rownames(fData(cds))[fData(cds)$use_for_ordering];de<-differentialGeneTest(cds[genes,],fullModelFormulaStr='~sm.ns(Pseudotime,df=3)',reducedModelFormulaStr='~1',cores=4);write.table(de,file.path(od,'pseudotime_dynamic_genes.tsv'),sep='\t',quote=FALSE,col.names=NA);saveRDS(cds,file.path(od,'ordered_monocle.rds'));writeLines('completed',file.path(od,'COMPLETE.txt'))
 },error=function(e){writeLines(c(conditionMessage(e),capture.output(sys.calls())),file.path(od,'error.txt'));message(branch,': ',conditionMessage(e))})}
writeLines(capture.output(sessionInfo()),'OUTPUT/sessionInfo.txt')
'''.replace('LIB',str(P/'R_library')).replace('OLD',str(D)).replace('OUTPUT',str(O))
(O/'run.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
with open(O/'R.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run.R')],stdout=log,stderr=subprocess.STDOUT,env=env)
v={'returncode':z.returncode,'branches':{b:{'completed':(O/b/'COMPLETE.txt').exists(),'error':(O/b/'error.txt').read_text() if (O/b/'error.txt').exists() else ''} for b in ['F3','F5']},'root_rule':'Original CytoTRACE-assisted Monocle2; reconstruct root as highest equal-donor mean CytoTRACE state; exploratory orientation, not fibrosis-to-HCC lineage proof. F3 mixedmesenchymal, not confirmed pureHSC; F3 conflictD9D10 excluded and recomputed.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
