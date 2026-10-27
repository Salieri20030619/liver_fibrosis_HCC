from pathlib import Path
import subprocess,os,json,time,pandas as pd
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';I=P/'continuation_followup_inputs_v1/F5';O=P/'continuation_expanded_irGSEA_v1';O.mkdir(exist_ok=False)
for i in range(180):
 if (I.parent/'validation.json').exists():break
 time.sleep(10)
else:raise RuntimeError('Inputs not ready')
libs=[P/'R_library',D/'irGSEA_dependencies_v1/R_project_library_v1',D/'trajectory_dependencies_v1/R_project_library_v1'];env=os.environ.copy();env['R_LIBS']=':'.join(map(str,libs));n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
r=r'''library(irGSEA);library(Seurat);library(Matrix);set.seed(20261020)
x<-as(readMM('INPUT/raw_counts.mtx'),'CsparseMatrix');rownames(x)<-readLines('INPUT/genes.txt');colnames(x)<-readLines('INPUT/cells.txt');m<-read.delim('INPUT/metadata.tsv',row.names=1,check.names=FALSE);stopifnot(identical(colnames(x),rownames(m)));x<-x%*%Diagonal(x=1e4/colSums(x));colnames(x)<-rownames(m);x@x<-log1p(x@x)
gmt<-strsplit(readLines('CURRENT/Hallmark_2025.1_Hs.gmt'),'\t');sets<-lapply(gmt,function(z)unique(z[-c(1,2)]));names(sets)<-sapply(gmt,function(z)z[1]);meth<-c('AUCell','UCell','singscore','ssgsea');writeLines('four_method_cell_scoring','OUTPUT/stage.txt')
obj<-irGSEA::irGSEA.score(object=x,min.cells=0,min.feature=0,seeds=20261020,ncores=4,custom=TRUE,geneset=sets,msigdb=FALSE,method=meth,kcdf='Gaussian',minGSSize=1,maxGSSize=Inf);saveRDS(obj,'OUTPUT/irGSEA_cell_scores.rds');stopifnot(all(meth %in% Assays(obj)));writeLines('record_state_aggregation_and_integrate','OUTPUT/stage.txt');records<-list();nr<-0;tests<-list();nt<-0
for(state in c('all_myeloid',sort(unique(m$refined_cluster)))){
 mm<-if(state=='all_myeloid')m else m[m$refined_cluster==state,,drop=FALSE];counts<-table(mm$donor);ids<-names(counts)[counts>=10];mm<-mm[mm$donor %in% ids,,drop=FALSE];if(length(ids)<6)next
 dm<-mm[match(ids,mm$donor),c('donor','study','F'),drop=FALSE];rownames(dm)<-ids;dm$group<-ifelse(dm$F>=3,'F3_F4','F0_F1_F2');avgs<-list()
 for(method in meth){scores<-GetAssayData(obj,assay=method,layer='scale.data');avg<-sapply(ids,function(id)rowMeans(scores[,rownames(mm)[mm$donor==id],drop=FALSE]));rownames(avg)<-rownames(scores);colnames(avg)<-ids;stopifnot(all(is.finite(avg)));avgs[[method]]<-avg;write.table(avg,paste0('OUTPUT/',state,'_',method,'_record_scores.tsv'),sep='\t',quote=FALSE,col.names=NA)}
 write.table(dm,paste0('OUTPUT/',state,'_record_metadata.tsv'),sep='\t',quote=FALSE,col.names=NA)
 for(study in c('pooled',sort(unique(dm$study)))){
  md<-if(study=='pooled')dm else dm[dm$study==study,,drop=FALSE];nc<-sum(md$group=='F3_F4');nn<-sum(md$group=='F0_F1_F2');key<-paste(state,study,sep='__');od<-file.path('OUTPUT',key);dir.create(od);nr<-nr+1
  if(min(nc,nn)<3){records[[nr]]<-data.frame(state=state,study=study,high_records=nc,low_records=nn,status='not_estimable_fewer3_pergroup',error='');next}
  tiny<-matrix(1,nrow=2,ncol=nrow(md),dimnames=list(c('placeholder1','placeholder2'),rownames(md)));don<-CreateSeuratObject(counts=tiny,meta.data=md,min.cells=0,min.features=0)
  for(method in meth){avg<-avgs[[method]][,rownames(md),drop=FALSE];ass<-CreateAssayObject(data=avg);ass<-SetAssayData(ass,layer='scale.data',new.data=avg);ass<-AddMetaData(ass,metadata=obj[[method]][[]]);don[[method]]<-ass
   pv<-apply(avg,1,function(v)wilcox.test(v[md$group=='F3_F4'],v[md$group=='F0_F1_F2'],exact=FALSE)$p.value);pv[!is.finite(pv)]<-1;nt<-nt+1;tests[[nt]]<-data.frame(state=state,study=study,method=method,pathway=rownames(avg),effect=rowMeans(avg[,md$group=='F3_F4',drop=FALSE])-rowMeans(avg[,md$group=='F0_F1_F2',drop=FALSE]),p=pv,BH_q=p.adjust(pv,'BH'),high_records=nc,low_records=nn)
  }
  saveRDS(don,file.path(od,'record_score_object.rds'));err<-'';res<-tryCatch(irGSEA::irGSEA.integrate(don,group.by='group',method=meth,p.value='p_val_adj'),error=function(e){err<<-conditionMessage(e);NULL});status<-if(is.null(res))'integration_no_result_requires_error_classification' else 'integrated';if(!is.null(res)){saveRDS(res,file.path(od,'integrated.rds'));writeLines(capture.output(str(res,max.level=2)),file.path(od,'integration_structure.txt'));for(nm in names(res)){if(is.data.frame(res[[nm]]))write.table(res[[nm]],file.path(od,paste0(nm,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)}}else writeLines(err,file.path(od,'integration_error.txt'));records[[nr]]<-data.frame(state=state,study=study,high_records=nc,low_records=nn,status=status,error=err);write.table(do.call(rbind,records),'OUTPUT/integration_status.tsv',sep='\t',quote=FALSE,row.names=FALSE)
 }
}
if(length(tests))write.table(do.call(rbind,tests),'OUTPUT/record_level_pathway_tests.tsv',sep='\t',quote=FALSE,row.names=FALSE);write.table(do.call(rbind,records),'OUTPUT/integration_status.tsv',sep='\t',quote=FALSE,row.names=FALSE);writeLines(capture.output(sessionInfo()),'OUTPUT/sessionInfo.txt');writeLines('completed','OUTPUT/COMPLETE.txt');writeLines('completed','OUTPUT/stage.txt')
'''.replace('INPUT',str(I)).replace('OUTPUT',str(O)).replace('CURRENT',str(P));(O/'run.R').write_text(r)
with open(O/'R.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run.R')],stdout=log,stderr=subprocess.STDOUT,env=env)
v={'status':'completed_scoring_tests_and_integration_attempts' if (O/'COMPLETE.txt').exists() else 'failed','R_returncode':z.returncode,'method':'Official irGSEA four methods AUCell/UCell/singscore/ssgsea, cached50Hallmark. All5010supportedmyeloidcells; record-state mean scores with >=10cells; perstudy and pooled F3/F4 versusF0/1/2 if>=3records/group.','limits':'Not undisclosedcustomfunction/metabolismsets. Source records not confirmedindependent. Empty significantsets must be distinguished from actual runtime errors. No cell-level pseudoreplication in stage contrasts.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
