library(monocle);set.seed(20261020)
sources<-c(F5='/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/end_to_end_20261020/continuation_expanded_trajectories_v1/F5',F3='/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/end_to_end_20261020/continuation_F3_CytoTRACE_serial_resume_v2')
rows<-list();k<-0
for(branch in names(sources)){
 src<-sources[[branch]];cds<-readRDS(file.path(src,'dispersion_checkpoint.rds'));old<-read.delim(file.path(src,'pseudotime_dynamic_genes.tsv'),row.names=1,check.names=FALSE);failed<-rownames(old)[old$status=='FAIL']
 stopifnot(all(is.finite(pData(cds)$Pseudotime)),all(is.finite(pData(cds)$Size_Factor)),all(pData(cds)$Size_Factor>0))
 for(g in failed){
  k<-k+1;cat(branch,g,k,'\n');flush.console();warn<-character();err<-''
  printed<-capture.output(ans<-tryCatch(withCallingHandlers(differentialGeneTest(cds[g,],fullModelFormulaStr='~sm.ns(Pseudotime,df=3)',reducedModelFormulaStr='~1',cores=1,verbose=TRUE),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL}))
  writeLines(c(printed,warn,err),file.path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/end_to_end_20261020/nonml_trajectory_failures_20261022',paste0(branch,'_',g,'_diagnostic.txt')))
  xx<-as.numeric(exprs(cds)[g,]);rows[[k]]<-data.frame(branch=branch,gene=g,old_status='FAIL',retry_status=if(is.null(ans))'OUTER_ERROR' else as.character(ans$status[1]),retry_p=if(is.null(ans))NA_real_ else ans$pval[1],cells=length(xx),detected_cells=sum(xx>0),unique_raw_values=length(unique(xx)),max_count=max(xx),message=paste(c(printed,warn,err),collapse=' | '))
  write.table(do.call(rbind,rows),'/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/end_to_end_20261020/nonml_trajectory_failures_20261022/failed_gene_retry.tsv',sep='\t',quote=TRUE,row.names=FALSE)
 }
 rm(cds);gc()
}
writeLines(capture.output(sessionInfo()),'/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/end_to_end_20261020/nonml_trajectory_failures_20261022/sessionInfo.txt');writeLines('completed','/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/end_to_end_20261020/nonml_trajectory_failures_20261022/COMPLETE.txt')
