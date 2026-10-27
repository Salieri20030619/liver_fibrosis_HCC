args<-commandArgs(trailingOnly=TRUE);root<-normalizePath(args[1]);job<-args[2]
P<-dirname(root);D<-file.path(dirname(P),"analysis_stage1_20261017_075157UTC")
.libPaths(c(file.path(root,"R_library"),file.path(P,"R_library"),file.path(D,"trajectory_dependencies_v1/R_project_library_v1"),file.path(D,"irGSEA_dependencies_v1/R_project_library_v1"),.libPaths()))
options(mc.cores=1,warn=1)
source(file.path(root,"code","repaired_core.R"))
methods<-read.delim(file.path(root,"code","frozen_113_methods.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(methods)==113,!anyDuplicated(methods$configuration))
X<-as.matrix(read.delim(file.path(root,"inputs","train_features.tsv"),row.names=1,check.names=FALSE))
meta<-read.delim(file.path(root,"inputs","train_design.tsv"),row.names=1,check.names=FALSE)
stopifnot(identical(rownames(X),rownames(meta)),nrow(X)==206,ncol(X)==30,all(meta$y==as.integer(meta$F>=3)),all(meta$F %in% 0:4))
y<-meta$y;folds<-read.delim(file.path(root,"inputs","outer_cv_folds.tsv"),row.names=1)
stopifnot(identical(rownames(folds),rownames(X)))
variants<-make_variant_table(methods)
out<-file.path(root,"jobs",job);dir.create(out,recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
if(job=="smoke") {
 set.seed(20261020)
 tr<-unlist(lapply(split(seq_along(y),y),function(ii)head(ii,40)))
 te<-setdiff(seq_along(y),tr)[1:20]
 # One parameter for each direct classifier, and one representative from each selector.
 mids<-methods$method_id[methods$selector=="None"]
 mids<-c(mids,sapply(setdiff(unique(methods$selector),"None"),function(s)methods$method_id[methods$selector==s][1]))
 vv<-variants[variants$method_id %in% mids,];vv<-vv[!duplicated(vv$method_id),]
 # Interior regularization value; other methods use first predeclared grid value.
 for(i in which(vv$classifier %in% c("Lasso","Ridge")|startsWith(vv$classifier,"Enet["))) {
  cp<-classifier_params(vv$classifier[i])[8];vv$classifier_parameter[i]<-cp;vv$key[i]<-variant_key(vv$method_id[i],vv$selector_parameter[i],cp)
 }
 z<-run_partition(X[tr,,drop=FALSE],y[tr],X[te,,drop=FALSE],methods,20261020,file.path(out,"representative_algorithms"),vv$key)
 write.table(z$status,file.path(out,"representative_status.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
 # Metamorphic test: transformation fitting cannot see heldout values.
 prep<-preprocess_fit(X[tr,,drop=FALSE]);changed<-X[te,,drop=FALSE]+1000
 stopifnot(identical(prep,preprocess_fit(X[tr,,drop=FALSE])))
 stopifnot(max(abs(preprocess_apply(prep,X[te,,drop=FALSE])-preprocess_apply(prep,changed)))>1)
 # Deliberately reordered columns must fail instead of using positional fallback.
 err<-try(preprocess_apply(prep,X[te,rev(colnames(X)),drop=FALSE]),silent=TRUE);stopifnot(inherits(err,"try-error"))
 # Positive AUROC direction and binary label convention are fixed, including inverse scores.
 stopifnot(auc_fixed(c(0,0,1,1),c(.1,.2,.8,.9))==1,auc_fixed(c(0,0,1,1),c(.9,.8,.2,.1))==0)
 # Forward selection must be able to leave the null model.
 set.seed(789);a<-matrix(rnorm(240),80,3,dimnames=list(paste0("s",1:80),paste0("g",1:3)))
 sy<-as.integer(a[,1]+rnorm(80,sd=.1)>0);sf<-step_fit(a,sy,"forward");stopifnot(length(coef(sf))>1)
 # Verify SVM scores are actual continuous P(y=1) and agree exactly with native prediction.
 if(requireNamespace("e1071",quietly=TRUE)) {
  aa<-preprocess_apply(prep,X[tr,,drop=FALSE]);bb<-preprocess_apply(prep,X[te,,drop=FALSE])
  fit<-fit_family("SVM",aa,y[tr],123);pp<-predict_family_one(fit,bb,"default")
  np<-attr(predict(fit$model,bb,probability=TRUE),"probabilities")[,"1"]
  stopifnot(max(abs(pp-np))<1e-12,length(unique(pp))>2)
 }
 v<-list(status="smoke_executed",attempted=nrow(z$status),completed=sum(z$status$status=="completed"),
  failures=z$status[z$status$status!="completed",],metamorphic_preprocessing=TRUE,ID_order_rejection=TRUE,
  fixed_AUC_direction=TRUE,forward_scope_test=TRUE,external_data_not_loaded=TRUE)
 jsonlite::write_json(v,file.path(out,"validation.json"),auto_unbox=TRUE,pretty=TRUE,na="null")
 quit(status=0)
}
if(job=="full") {
 ids<-seq_len(nrow(X));inner<-data.frame(sample=rownames(X),inner_fold=folds$fold);seed<-20262020
} else {
 k<-as.integer(sub("outer_","",job));ids<-which(folds$fold!=k)
 inner<-read.delim(file.path(root,"inputs",paste0("inner_folds_",k,".tsv")))
 stopifnot(identical(inner$sample,rownames(X)[ids]));seed<-20261020+k*1000
}
xx<-X[ids,,drop=FALSE];yy<-y[ids]
pred<-matrix(NA_real_,nrow(xx),nrow(variants),dimnames=list(rownames(xx),variants$key))
for(k in sort(unique(inner$inner_fold))) {
 tr<-which(inner$inner_fold!=k);te<-which(inner$inner_fold==k);od<-file.path(out,paste0("inner_",k))
 if(file.exists(file.path(od,"DONE.txt"))) {
  pp<-readRDS(file.path(od,"predictions.rds"))
 } else {
  zz<-run_partition(xx[tr,,drop=FALSE],yy[tr],xx[te,,drop=FALSE],methods,seed+k,od)
  pp<-zz$predictions;writeLines("completed",file.path(od,"DONE.txt"))
 }
 stopifnot(identical(rownames(pp),rownames(xx)[te]),identical(colnames(pp),variants$key))
 pred[te,]<-pp
 cat("INNER_COMPLETE",job,k,"\n");flush.console()
}
saveRDS(pred,file.path(out,"inner_OOF_predictions.rds"))
scores<-variants;scores$all_folds_complete<-apply(pred,2,function(p)all(is.finite(p)))
scores$mean_fold_AUC<-NA_real_;scores$pooled_AUC<-NA_real_;scores$Brier<-NA_real_
for(j in which(scores$all_folds_complete)) {
 pp<-pred[,j];scores$mean_fold_AUC[j]<-mean(sapply(sort(unique(inner$inner_fold)),function(k)auc_fixed(yy[inner$inner_fold==k],pp[inner$inner_fold==k])))
 scores$pooled_AUC[j]<-auc_fixed(yy,pp);scores$Brier[j]<-mean((yy-pp)^2)
}
write.table(scores,file.path(out,"all_variant_inner_metrics.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
best<-list()
for(id in methods$method_id) {
 z<-scores[scores$method_id==id & scores$all_folds_complete,,drop=FALSE]
 if(!nrow(z))next
 z<-z[order(-z$mean_fold_AUC,z$Brier,match(z$key,variants$key)),,drop=FALSE];best[[id]]<-z[1,,drop=FALSE]
}
if(!length(best))stop("No configuration completed all inner folds")
best<-do.call(rbind,best);best<-best[order(-best$mean_fold_AUC,best$Brier,match(best$method_id,methods$method_id)),,drop=FALSE]
write.table(best,file.path(out,"best_parameters_per_configuration.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
chosen<-best[1,,drop=FALSE];jsonlite::write_json(chosen,file.path(out,"TRAINING_ONLY_SELECTION_FROZEN.json"),dataframe="rows",auto_unbox=TRUE,pretty=TRUE)
if(job=="full") {
 # The external matrix is deliberately first opened AFTER the training-only decision is persisted.
 ee<-as.matrix(read.delim(file.path(root,"inputs","external_features.tsv"),row.names=1,check.names=FALSE))
 em<-read.delim(file.path(root,"inputs","external_design.tsv"),row.names=1,check.names=FALSE)
 stopifnot(identical(rownames(ee),rownames(em)),identical(colnames(ee),colnames(X)),nrow(ee)==112,all(em$y==as.integer(em$F>=3)))
 z<-run_partition(X,y,ee,methods,seed+900,file.path(out,"frozen_external"),requested_keys=best$key,save_key=chosen$key)
 p<-pred[,chosen$key];thresholds<-sort(unique(c(0,p,1)));youden<-sapply(thresholds,function(t)mean(p[yy==1]>=t)+mean(p[yy==0]<t)-1)
 threshold<-thresholds[which.max(youden)]
 write.table(data.frame(sample=rownames(X),fold=inner$inner_fold,y=y,probability=p),file.path(out,"selected_training_OOF.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 write.table(data.frame(sample=rownames(ee),y=em$y,z$predictions,check.names=FALSE),gzfile(file.path(out,"all_frozen_external_predictions.tsv.gz")),sep="\t",quote=FALSE,row.names=FALSE)
 jsonlite::write_json(list(status="completed",selected_key=chosen$key,selected_method=chosen$configuration,
  threshold_from_training_OOF=threshold,eligible_configurations=nrow(best),external_completed=sum(z$status$status=="completed"),
  no_external_selection=TRUE),file.path(out,"validation.json"),auto_unbox=TRUE,pretty=TRUE)
} else {
 hold<-which(folds$fold==as.integer(sub("outer_","",job)))
 z<-run_partition(xx,yy,X[hold,,drop=FALSE],methods,seed+900,file.path(out,"outer_prediction"),requested_keys=best$key)
 write.table(data.frame(sample=rownames(X)[hold],y=y[hold],z$predictions,check.names=FALSE),gzfile(file.path(out,"all_configuration_outer_predictions.tsv.gz")),sep="\t",quote=FALSE,row.names=FALSE)
 jsonlite::write_json(list(status="completed",outer_fold=as.integer(sub("outer_","",job)),selected_key=chosen$key,
  selected_method=chosen$configuration,eligible_configurations=nrow(best),outer_completed=sum(z$status$status=="completed")),
  file.path(out,"validation.json"),auto_unbox=TRUE,pretty=TRUE)
}
writeLines("completed",file.path(out,"COMPLETE.txt"))
