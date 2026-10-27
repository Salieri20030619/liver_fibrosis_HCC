from pathlib import Path
import subprocess,os,json,time
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'nonml_trajectory_failures_20261022';O.mkdir(exist_ok=False)
env=os.environ.copy();env['R_LIBS']=':'.join(map(str,[P/'R_monocle_legacy',P/'R_library',D/'trajectory_dependencies_v1/R_project_library_v1',D/'irGSEA_dependencies_v1/R_project_library_v1']))
n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
r=r'''library(monocle);set.seed(20261020)
sources<-c(F5='ROOT/continuation_expanded_trajectories_v1/F5',F3='ROOT/continuation_F3_CytoTRACE_serial_resume_v2')
rows<-list();k<-0
for(branch in names(sources)){
 src<-sources[[branch]];cds<-readRDS(file.path(src,'dispersion_checkpoint.rds'));old<-read.delim(file.path(src,'pseudotime_dynamic_genes.tsv'),row.names=1,check.names=FALSE);failed<-rownames(old)[old$status=='FAIL']
 stopifnot(all(is.finite(pData(cds)$Pseudotime)),all(is.finite(pData(cds)$Size_Factor)),all(pData(cds)$Size_Factor>0))
 for(g in failed){
  k<-k+1;cat(branch,g,k,'\n');flush.console();warn<-character();err<-''
  printed<-capture.output(ans<-tryCatch(withCallingHandlers(differentialGeneTest(cds[g,],fullModelFormulaStr='~sm.ns(Pseudotime,df=3)',reducedModelFormulaStr='~1',cores=1,verbose=TRUE),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);NULL}))
  writeLines(c(printed,warn,err),file.path('OUTPUT',paste0(branch,'_',g,'_diagnostic.txt')))
  xx<-as.numeric(exprs(cds)[g,]);rows[[k]]<-data.frame(branch=branch,gene=g,old_status='FAIL',retry_status=if(is.null(ans))'OUTER_ERROR' else as.character(ans$status[1]),retry_p=if(is.null(ans))NA_real_ else ans$pval[1],cells=length(xx),detected_cells=sum(xx>0),unique_raw_values=length(unique(xx)),max_count=max(xx),message=paste(c(printed,warn,err),collapse=' | '))
  write.table(do.call(rbind,rows),'OUTPUT/failed_gene_retry.tsv',sep='\t',quote=TRUE,row.names=FALSE)
 }
 rm(cds);gc()
}
writeLines(capture.output(sessionInfo()),'OUTPUT/sessionInfo.txt');writeLines('completed','OUTPUT/COMPLETE.txt')
'''.replace('ROOT',str(P)).replace('OUTPUT',str(O));(O/'run.R').write_text(r)
start=time.time()
with open(O/'R.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run.R')],env=env,stdout=log,stderr=subprocess.STDOUT)
(O/'validation.json').write_text(json.dumps({'returncode':z.returncode,'seconds':time.time()-start,'completed':(O/'COMPLETE.txt').exists(),'scope':'Only original failed genes; identical model and cached trajectory/dispersion; serial verbose retry. No original result overwritten; retry q-values not used as full-family BH.'},indent=2))
