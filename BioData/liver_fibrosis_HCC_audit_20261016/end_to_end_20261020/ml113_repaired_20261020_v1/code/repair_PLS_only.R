args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);P<-dirname(root);D<-file.path(dirname(P),"analysis_stage1_20261017_075157UTC")
.libPaths(c(file.path(root,"R_randomForestSRC_3.3.2"),file.path(root,"R_library"),file.path(P,"R_library"),file.path(D,"trajectory_dependencies_v1/R_project_library_v1"),file.path(D,"irGSEA_dependencies_v1/R_project_library_v1"),.libPaths()))
suppressPackageStartupMessages(library(plsRglm));options(warn=1)
source(file.path(root,"code","repaired_core_PLS_v5.R"))
methods<-read.delim(file.path(root,"code","frozen_113_methods.tsv"),check.names=FALSE)
variants<-make_variant_table(methods);plsids<-methods$method_id[methods$classifier=="plsRglm"];stopifnot(length(plsids)==6)
plkeys<-variants$key[variants$classifier=="plsRglm"];otherkeys<-setdiff(variants$key,plkeys)
X<-as.matrix(read.delim(file.path(root,"inputs","train_features.tsv"),row.names=1,check.names=FALSE));m<-read.delim(file.path(root,"inputs","train_design.tsv"),row.names=1);y<-m$y
fold<-read.delim(file.path(root,"inputs","outer_cv_folds.tsv"),row.names=1)$fold
stopifnot(identical(rownames(X),rownames(m)))
base<-file.path(root,"PLS_multicomponent_repair_v2");dir.create(base);dir.create(file.path(base,"jobs"))
# Test BOTH naming conventions with the same sparse PLS algorithm.
tr<-unlist(lapply(split(seq_along(y),y),function(ii)head(ii,40)));te<-setdiff(seq_along(y),tr)[1:20]
prep<-preprocess_fit(X[tr,]);aa<-preprocess_apply(prep,X[tr,]);bb<-preprocess_apply(prep,X[te,])
test<-list()
for(nt in c("1","2")) {
 f<-fit_family("plsRglm",aa,y[tr],123,nt);obj<-f$model[[nt]];stopifnot(!inherits(obj,"fit_error"),obj$computed_nt==as.integer(nt))
 scores<-predict(obj,newdata=as.data.frame(aa),type="scores",verbose=FALSE)
 projection_error<-max(abs(unname(as.matrix(scores))-unname(as.matrix(obj$tt))))
 p<-predict_family_one(f,aa,nt);identity_error<-max(abs(p-as.numeric(obj$FinalModel$fitted.values)))
 stopifnot(projection_error<1e-8,identity_error<1e-8,all(is.finite(predict_family_one(f,bb,nt))))
 test[[nt]]<-list(components=as.integer(nt),projection_error=projection_error,training_fitted_probability_error=identity_error)
}
jsonlite::write_json(test,file.path(base,"multi_component_acceptance.json"),auto_unbox=TRUE,pretty=TRUE)
for(job in c(paste0("outer_",0:9),"full")) {
 old<-file.path(root,"jobs",job);out<-file.path(base,"jobs",job);dir.create(out)
 stopifnot(file.exists(file.path(old,"COMPLETE.txt")))
 if(job=="full"){ids<-seq_len(nrow(X));inner<-data.frame(sample=rownames(X),inner_fold=fold)}
 else {k<-as.integer(sub("outer_","",job));ids<-which(fold!=k);inner<-read.delim(file.path(root,"inputs",paste0("inner_folds_",k,".tsv")))}
 xx<-X[ids,,drop=FALSE];yy<-y[ids];stopifnot(identical(rownames(xx),inner$sample))
 previous<-readRDS(file.path(old,"inner_OOF_predictions.rds"));pred<-previous
 for(k in sort(unique(inner$inner_fold))) {
  oi<-file.path(old,paste0("inner_",k));ni<-file.path(out,paste0("inner_",k));dir.create(ni)
  pp<-readRDS(file.path(oi,"predictions.rds"));before<-pp
  st<-read.delim(file.path(oi,"variant_status.tsv"),check.names=FALSE,na.strings="");st$error[is.na(st$error)]<-""
  ii<-read.delim(file.path(oi,"partition_IDs.tsv"));trainid<-ii$sample[ii$role=="train"];testid<-ii$sample[ii$role=="heldout"]
  prep<-readRDS(file.path(oi,"preprocessing.rds"));stopifnot(identical(prep$training_ids,trainid))
  a<-preprocess_apply(prep,X[trainid,,drop=FALSE]);b<-preprocess_apply(prep,X[testid,,drop=FALSE]);yt<-y[match(trainid,rownames(X))]
  pk<-which(st$classifier=="plsRglm")
  for(s in unique(st$selector[pk]))for(sp in unique(st$selector_parameter[pk][st$selector[pk]==s])) {
   kk<-pk[st$selector[pk]==s & st$selector_parameter[pk]==sp]
   if(all(st$status[kk]=="selector_failed"))next
   fs<-strsplit(st$features[kk[1]],";",fixed=TRUE)[[1]];stopifnot(length(fs)>0,all(fs %in% colnames(a)))
   fit<-fit_family("plsRglm",a[,fs,drop=FALSE],yt,20261020,st$classifier_parameter[kk])
   for(j in kk) {
    z<-tryCatch(predict_family_one(fit,b[,fs,drop=FALSE],st$classifier_parameter[j]),error=function(e)e)
    if(inherits(z,"error")) {pp[,st$key[j]]<-NA_real_;st$status[j]<-"prediction_failed";st$error[j]<-conditionMessage(z)}
    else {pp[,st$key[j]]<-z;st$status[j]<-"completed";st$error[j]<-""}
   }
  }
  stopifnot(identical(before[,otherkeys,drop=FALSE],pp[,otherkeys,drop=FALSE]))
  pred[match(testid,rownames(pred)),]<-pp
  saveRDS(pp,file.path(ni,"predictions.rds"));write.table(st,file.path(ni,"variant_status.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  file.copy(file.path(oi,"partition_IDs.tsv"),file.path(ni,"partition_IDs.tsv"));file.copy(file.path(oi,"preprocessing.rds"),file.path(ni,"preprocessing.rds"))
  writeLines("completed",file.path(ni,"DONE.txt"))
 }
 stopifnot(identical(previous[,otherkeys,drop=FALSE],pred[,otherkeys,drop=FALSE]))
 saveRDS(pred,file.path(out,"inner_OOF_predictions.rds"))
 scores<-variants;scores$all_folds_complete<-apply(pred,2,function(p)all(is.finite(p)));scores$mean_fold_AUC<-NA_real_;scores$pooled_AUC<-NA_real_;scores$Brier<-NA_real_
 for(j in which(scores$all_folds_complete)) {
  p<-pred[,j];scores$mean_fold_AUC[j]<-mean(sapply(sort(unique(inner$inner_fold)),function(k)auc_fixed(yy[inner$inner_fold==k],p[inner$inner_fold==k])));scores$pooled_AUC[j]<-auc_fixed(yy,p);scores$Brier[j]<-mean((yy-p)^2)
 }
 write.table(scores,file.path(out,"all_variant_inner_metrics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 best<-list()
 for(id in methods$method_id) {
  z<-scores[scores$method_id==id & scores$all_folds_complete,,drop=FALSE]
  if(!nrow(z))next
  z<-z[order(-z$mean_fold_AUC,z$Brier,match(z$key,variants$key)),,drop=FALSE];best[[id]]<-z[1,,drop=FALSE]
 }
 best<-do.call(rbind,best);best<-best[order(-best$mean_fold_AUC,best$Brier,match(best$method_id,methods$method_id)),,drop=FALSE];chosen<-best[1,,drop=FALSE]
 priorbest<-read.delim(file.path(old,"best_parameters_per_configuration.tsv"),check.names=FALSE)
 for(id in setdiff(methods$method_id,plsids))stopifnot(identical(best$key[best$method_id==id],priorbest$key[priorbest$method_id==id]))
 write.table(best,file.path(out,"best_parameters_per_configuration.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 jsonlite::write_json(chosen,file.path(out,"TRAINING_ONLY_SELECTION_FROZEN.json"),dataframe="rows",auto_unbox=TRUE,pretty=TRUE)
 od<-if(job=="full")"frozen_external" else "outer_prediction";op<-file.path(old,od);np<-file.path(out,od);dir.create(np)
 # All training-only selections have been re-frozen before reading external predictions/data.
 oldp<-readRDS(file.path(op,"predictions.rds"));oldst<-read.delim(file.path(op,"variant_status.tsv"),check.names=FALSE,na.strings="")
 prep<-readRDS(file.path(op,"preprocessing.rds"));stopifnot(identical(prep$training_ids,rownames(xx)))
 if(job=="full") {
  rawb<-as.matrix(read.delim(file.path(root,"inputs","external_features.tsv"),row.names=1,check.names=FALSE))
  ey<-read.delim(file.path(root,"inputs","external_design.tsv"),row.names=1)$y
 } else rawb<-X[fold==as.integer(sub("outer_","",job)),,drop=FALSE]
 a<-preprocess_apply(prep,xx);b<-preprocess_apply(prep,rawb);stopifnot(identical(rownames(b),rownames(oldp)))
 pp<-matrix(NA_real_,nrow(b),nrow(best),dimnames=list(rownames(b),best$key));statuses<-list()
 for(i in seq_len(nrow(best))) {
  row<-best[i,];key<-row$key
  if(row$classifier!="plsRglm") {
   stopifnot(key %in% colnames(oldp));pp[,key]<-oldp[,key];statuses[[key]]<-oldst[oldst$key==key,,drop=FALSE];next
  }
  z<-oldst[oldst$selector==row$selector & oldst$selector_parameter==row$selector_parameter,,drop=FALSE]
  # Reuse existing selector outputs. Only deterministic non-RF selector parameters
  # absent from the earlier chosen set may need fitting on this training partition.
  fs<-tryCatch({
   if(nrow(z)&&!is.na(z$features[1])&&nzchar(z$features[1]))strsplit(z$features[1],";",fixed=TRUE)[[1]]
   else {stopifnot(row$selector!="RF");select_features(row$selector,row$selector_parameter,a,yy,20261020)}
  },error=function(e)e)
  sr<-variants[variants$key==key,];sr$status<-"not_run";sr$error<-"";sr$feature_count<-NA_integer_;sr$features<-""
  if(inherits(fs,"error")){sr$status<-"selector_failed";sr$error<-conditionMessage(fs)}
  else {
   sr$feature_count<-length(fs);sr$features<-paste(fs,collapse=";")
   fit<-fit_family("plsRglm",a[,fs,drop=FALSE],yy,20261020,row$classifier_parameter)
   predone<-tryCatch(predict_family_one(fit,b[,fs,drop=FALSE],row$classifier_parameter),error=function(e)e)
   if(inherits(predone,"error")){sr$status<-"prediction_failed";sr$error<-conditionMessage(predone)}
   else {
    pp[,key]<-predone;sr$status<-"completed"
    if(job=="full" && key==chosen$key)saveRDS(list(preprocessing=prep,selected_features=fs,model=fit,parameter=row$classifier_parameter,variant=row),file.path(np,"selected_frozen_model.rds"))
   }
  }
  statuses[[key]]<-sr
 }
 for(id in setdiff(methods$method_id,plsids)) {
  key<-best$key[best$method_id==id]
  if(length(key))stopifnot(identical(as.numeric(pp[,key]),as.numeric(oldp[,key])))
 }
 st<-do.call(rbind,statuses);write.table(st,file.path(np,"variant_status.tsv"),sep="\t",quote=FALSE,row.names=FALSE);saveRDS(pp,file.path(np,"predictions.rds"));saveRDS(prep,file.path(np,"preprocessing.rds"));file.copy(file.path(op,"partition_IDs.tsv"),file.path(np,"partition_IDs.tsv"))
 if(job=="full") {
  if(chosen$classifier!="plsRglm") {
   prev<-jsonlite::read_json(file.path(old,"validation.json"),simplifyVector=TRUE);stopifnot(prev$selected_key==chosen$key)
   file.copy(file.path(op,"selected_frozen_model.rds"),file.path(np,"selected_frozen_model.rds"))
  }
  p<-pred[,chosen$key];thresholds<-sort(unique(c(0,p,1)));j<-sapply(thresholds,function(t)mean(p[yy==1]>=t)+mean(p[yy==0]<t)-1);threshold<-thresholds[which.max(j)]
  write.table(data.frame(sample=rownames(X),fold=inner$inner_fold,y=y,probability=p),file.path(out,"selected_training_OOF.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(data.frame(sample=rownames(rawb),y=ey,pp,check.names=FALSE),gzfile(file.path(out,"all_frozen_external_predictions.tsv.gz")),sep="\t",quote=FALSE,row.names=FALSE)
  v<-list(status="completed",selected_key=chosen$key,selected_method=chosen$configuration,threshold_from_training_OOF=threshold,eligible_configurations=nrow(best),external_completed=sum(st$status=="completed"),no_external_selection=TRUE)
 } else {
  write.table(data.frame(sample=rownames(rawb),y=y[match(rownames(rawb),rownames(X))],pp,check.names=FALSE),gzfile(file.path(out,"all_configuration_outer_predictions.tsv.gz")),sep="\t",quote=FALSE,row.names=FALSE)
  v<-list(status="completed",outer_fold=as.integer(sub("outer_","",job)),selected_key=chosen$key,selected_method=chosen$configuration,eligible_configurations=nrow(best),outer_completed=sum(st$status=="completed"))
 }
 v$nonPLS107_predictions_reused_identically<-TRUE
 jsonlite::write_json(v,file.path(out,"validation.json"),auto_unbox=TRUE,pretty=TRUE);writeLines("completed",file.path(out,"COMPLETE.txt"))
 cat("REPAIRED",job,"\n");flush.console()
}
jsonlite::write_json(list(status="completed",repaired_configurations=6,reused_configurations=107,all_nonPLS_predictions_identical=TRUE,reason="PLS1 uses Comp_k terms; multicomponent uses tt.Comp_k. Explicit canonical name mapping with fitted-probability identity tests. No changed grids or external selection."),file.path(base,"validation.json"),auto_unbox=TRUE,pretty=TRUE)
