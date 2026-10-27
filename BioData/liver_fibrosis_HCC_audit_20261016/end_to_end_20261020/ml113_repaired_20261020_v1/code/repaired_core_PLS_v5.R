# Frozen 113 feature-selection -> classifier implementation.
# No external labels enter these functions. All preprocessing is fitted on a.
auc_fixed <- function(y,p) {
 stopifnot(length(y)==length(p),all(y %in% c(0,1)),length(unique(y))==2,
           all(is.finite(p)),all(p>=0 & p<=1))
 n1<-sum(y==1);n0<-sum(y==0)
 (sum(rank(p,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}
preprocess_fit <- function(a) {
 a<-as.matrix(a);a[!is.finite(a)]<-NA_real_
 med<-apply(a,2,median,na.rm=TRUE);if(any(!is.finite(med)))stop("All-missing training feature")
 for(j in seq_len(ncol(a)))a[is.na(a[,j]),j]<-med[j]
 mu<-colMeans(a);ss<-apply(a,2,sd);keep<-names(ss)[is.finite(ss)&ss>0]
 if(!length(keep))stop("No nonconstant training features")
 list(median=med,mean=mu,sd=ss,keep=keep,training_ids=rownames(a))
}
preprocess_apply <- function(z,a) {
 stopifnot(identical(colnames(a),names(z$median)))
 a<-as.matrix(a);a[!is.finite(a)]<-NA_real_
 for(j in seq_len(ncol(a)))a[is.na(a[,j]),j]<-z$median[j]
 scale(a[,z$keep,drop=FALSE],center=z$mean[z$keep],scale=z$sd[z$keep])
}
step_fit <- function(a,y,direction) {
 dat<-data.frame(y=y,a,check.names=FALSE)
 upper<-reformulate(colnames(a),response="y")
 if(direction=="forward") {
  fit<-stats::glm(y~1,data=dat,family=binomial())
  fit<-stats::step(fit,scope=list(lower=y~1,upper=upper),direction="forward",trace=0)
 } else {
  fit<-stats::glm(upper,data=dat,family=binomial())
  fit<-stats::step(fit,scope=list(lower=y~1,upper=upper),direction=direction,trace=0)
 }
 if(!isTRUE(fit$converged))stop("Final Stepglm fit did not converge")
 fit
}
selector_params <- function(s) switch(s,"None"="default","Lasso"=c("0.01","0.03","0.1"),
 "glmBoost"=c("40","70","100"),"RF"="default_var.select",
 "Stepglm[both]"="AIC","Stepglm[backward]"="AIC",stop("Unknown selector"))
select_features <- function(s,param,a,y,seed) {
 set.seed(seed)
 if(s=="None")return(colnames(a))
 if(s=="Lasso") {
  fit<-glmnet::glmnet(a,y,family="binomial",alpha=1,lambda=as.numeric(param),standardize=FALSE)
  cc<-as.matrix(coef(fit));g<-rownames(cc)[cc[,1]!=0]
 } else if(s=="glmBoost") {
  dat<-data.frame(y=factor(y,levels=c(0,1)),a,check.names=FALSE)
  fit<-mboost::glmboost(y~.,data=dat,family=mboost::Binomial(),control=mboost::boost_control(mstop=as.integer(param)))
  cc<-coef(fit);g<-names(cc)[abs(cc)>0]
 } else if(s=="RF") {
  fit<-randomForestSRC::rfsrc(y~.,data=data.frame(y=factor(y,levels=c(0,1)),a,check.names=FALSE),
    ntree=1000,nodesize=5,importance=TRUE,proximity=FALSE,forest=TRUE,nthread=1)
  g<-randomForestSRC::var.select(fit,verbose=FALSE)$topvars
 } else {
  direction<-if(s=="Stepglm[both]")"both" else "backward"
  fit<-step_fit(a,y,direction);g<-names(coef(fit))
 }
 g<-setdiff(g,c("(Intercept)","Intercept"))
 if(!length(g))stop("Empty selected feature set; no fallback")
 if(!all(g %in% colnames(a)) || anyDuplicated(g))stop("Invalid selected feature IDs")
 colnames(a)[colnames(a) %in% g]
}
classifier_params <- function(cl) {
 if(cl %in% c("Lasso","Ridge") || startsWith(cl,"Enet["))return(format(10^seq(0,-4,length.out=16),scientific=TRUE,digits=12,trim=TRUE))
 switch(cl,"glmBoost"=c("40","70","100"),"plsRglm"=as.character(1:10),
  "GBM"=c("500","2000","10000"),"XGBoost"=as.character(1:10),"default")
}
fit_family <- function(cl,a,y,seed,requested=NULL) {
 set.seed(seed);params<-classifier_params(cl)
 if(!is.null(requested))params<-params[params %in% requested]
 if(!length(params))stop("No requested parameter")
 dat<-as.data.frame(a,check.names=FALSE)
 if(cl %in% c("Lasso","Ridge") || startsWith(cl,"Enet[")) {
  alpha<-if(cl=="Lasso")1 else if(cl=="Ridge")0 else as.numeric(sub("Enet\\[alpha=(.*)\\]","\\1",cl))
  fit<-glmnet::glmnet(a,y,family="binomial",alpha=alpha,lambda=sort(unique(as.numeric(params)),decreasing=TRUE),standardize=FALSE)
 } else if(startsWith(cl,"Stepglm[")) {
  fit<-step_fit(a,y,sub("Stepglm\\[(.*)\\]","\\1",cl))
 } else if(cl=="SVM") {
  fit<-e1071::svm(x=a,y=factor(y,levels=c(0,1)),kernel="radial",cost=1,gamma=1/ncol(a),scale=FALSE,probability=TRUE)
 } else if(cl=="LDA") {
  fit<-MASS::lda(x=a,grouping=factor(y,levels=c(0,1)))
 } else if(cl=="glmBoost") {
  dat$y<-factor(y,levels=c(0,1));fit<-mboost::glmboost(y~.,data=dat,family=mboost::Binomial(),control=mboost::boost_control(mstop=max(as.integer(params))))
 } else if(cl=="plsRglm") {
  # Each valid nt is passed into the actual model; no unused CV result.
  fit<-setNames(lapply(params,function(p)tryCatch({
   nt<-as.integer(p);if(nt>min(ncol(a),nrow(a)-2))stop("Component count exceeds available rank bound")
   obj<-plsRglm::plsRglm(object=y,dataX=dat,nt=nt,modele="pls-glm-logistic",verbose=FALSE,sparse=TRUE)
   if(obj$computed_nt!=nt)stop("PLS extracted fewer components than requested")
   if(!isTRUE(obj$FinalModel$converged))stop("PLS final GLM did not converge")
   obj
  },error=function(e)structure(list(message=conditionMessage(e)),class="fit_error"))),params)
 } else if(cl=="RF") {
  dat$y<-factor(y,levels=c(0,1));fit<-randomForestSRC::rfsrc(y~.,data=dat,ntree=1000,nodesize=5,importance=TRUE,proximity=FALSE,forest=TRUE,nthread=1)
 } else if(cl=="GBM") {
  dat$y<-y;fit<-gbm::gbm(y~.,data=dat,distribution="bernoulli",n.trees=max(as.integer(params)),interaction.depth=3,
   n.minobsinnode=10,shrinkage=.001,bag.fraction=.5,cv.folds=0,n.cores=1,verbose=FALSE)
 } else if(cl=="XGBoost") {
  fit<-xgboost::xgb.train(params=list(max_depth=2,eta=1,nthread=1,objective="binary:logistic",seed=seed),
   data=xgboost::xgb.DMatrix(a,label=y),nrounds=max(as.integer(params)),verbose=0)
 } else if(cl=="NaiveBayes") {
  fit<-e1071::naiveBayes(x=a,y=factor(y,levels=c(0,1)))
 } else stop("Unknown classifier")
 list(classifier=cl,model=fit,params=params,features=colnames(a))
}
predict_family_one <- function(f,b,param) {
 stopifnot(identical(colnames(b),f$features));cl<-f$classifier;fit<-f$model;dat<-as.data.frame(b,check.names=FALSE)
 if(cl %in% c("Lasso","Ridge") || startsWith(cl,"Enet["))p<-predict(fit,newx=b,s=as.numeric(param),type="response")
 else if(startsWith(cl,"Stepglm["))p<-predict(fit,newdata=dat,type="response")
 else if(cl=="SVM") {
  z<-predict(fit,newdata=b,probability=TRUE);pp<-attr(z,"probabilities");stopifnot("1" %in% colnames(pp));p<-pp[,"1"]
 } else if(cl=="LDA")p<-predict(fit,newdata=b)$posterior[,"1"]
 else if(cl=="glmBoost")p<-predict(fit[as.integer(param)],newdata=dat,type="response")
 else if(cl=="plsRglm") {
  obj<-fit[[param]];if(inherits(obj,"fit_error"))stop(obj$message)
  # plsRglm1.7.1 response helper indexes Comp_ coefficients using different
  # generated design names and returns NA even on training data.
  # Keep the same frozen PLS projection and final GLM; use native GLM prediction.
  scores<-predict(obj,newdata=dat,type="scores",verbose=FALSE)
  scoredata<-as.data.frame(scores,check.names=FALSE)
  required<-all.vars(stats::delete.response(stats::terms(obj$FinalModel)))
  canonical<-sub("^tt\\.","",required)
  if(!setequal(colnames(scoredata),canonical)||anyDuplicated(canonical))stop("PLS score/model term name mismatch")
  scoredata<-scoredata[,canonical,drop=FALSE];colnames(scoredata)<-required
  p<-stats::predict(obj$FinalModel,newdata=scoredata,type="response")
 } else if(cl=="RF")p<-predict(fit,newdata=dat)$predicted[,"1"]
 else if(cl=="GBM")p<-predict(fit,newdata=dat,type="response",n.trees=as.integer(param))
 else if(cl=="XGBoost")p<-predict(fit,newdata=xgboost::xgb.DMatrix(b),iterationrange=c(1L,as.integer(param)))
 else if(cl=="NaiveBayes")p<-predict(fit,newdata=dat,type="raw")[,"1"]
 p<-as.numeric(p)
 if(length(p)!=nrow(b)||any(!is.finite(p))||any(p<0 | p>1))stop("Invalid continuous positive probabilities")
 p
}
variant_key <- function(id,sp,cp)paste(id,sp,cp,sep="|")
make_variant_table <- function(methods) {
 out<-list()
 for(i in seq_len(nrow(methods)))for(sp in selector_params(methods$selector[i]))for(cp in classifier_params(methods$classifier[i]))
  out[[length(out)+1]]<-data.frame(method_id=methods$method_id[i],configuration=methods$configuration[i],
   selector=methods$selector[i],classifier=methods$classifier[i],selector_parameter=sp,classifier_parameter=cp,
   key=variant_key(methods$method_id[i],sp,cp))
 do.call(rbind,out)
}
run_partition <- function(a,y,b,methods,seed,outdir,requested_keys=NULL,save_key=NULL) {
 stopifnot(nrow(a)==length(y),!anyDuplicated(rownames(a)),!anyDuplicated(rownames(b)),!length(intersect(rownames(a),rownames(b))))
 dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
 variants<-make_variant_table(methods)
 if(!is.null(requested_keys))variants<-variants[variants$key %in% requested_keys,,drop=FALSE]
 pr<-matrix(NA_real_,nrow(b),nrow(variants),dimnames=list(rownames(b),variants$key))
 status<-variants;status$status<-"not_run";status$error<-"";status$feature_count<-NA_integer_;status$features<-""
 prep<-preprocess_fit(a);aa<-preprocess_apply(prep,a);bb<-preprocess_apply(prep,b)
 saveRDS(prep,file.path(outdir,"preprocessing.rds"))
 write.table(data.frame(role=c(rep("train",nrow(a)),rep("heldout",nrow(b))),sample=c(rownames(a),rownames(b))),file.path(outdir,"partition_IDs.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 pairs<-unique(variants[,c("selector","selector_parameter")])
 for(j in seq_len(nrow(pairs))) {
  s<-pairs$selector[j];sp<-pairs$selector_parameter[j];ii<-which(variants$selector==s & variants$selector_parameter==sp)
  fs<-tryCatch(select_features(s,sp,aa,y,seed+j),error=function(e)e)
  if(inherits(fs,"error")) {status$status[ii]<-"selector_failed";status$error[ii]<-conditionMessage(fs);next}
  status$feature_count[ii]<-length(fs);status$features[ii]<-paste(fs,collapse=";")
  for(cl in unique(variants$classifier[ii])) {
   kk<-ii[variants$classifier[ii]==cl];cp<-unique(variants$classifier_parameter[kk])
   fit<-tryCatch(fit_family(cl,aa[,fs,drop=FALSE],y,seed+100+j,cp),error=function(e)e)
   if(inherits(fit,"error")){status$status[kk]<-"classifier_failed";status$error[kk]<-conditionMessage(fit);next}
   for(k in kk) {
    z<-tryCatch(predict_family_one(fit,bb[,fs,drop=FALSE],variants$classifier_parameter[k]),error=function(e)e)
    if(inherits(z,"error")){status$status[k]<-"prediction_failed";status$error[k]<-conditionMessage(z)}
    else {pr[,k]<-z;status$status[k]<-"completed"}
    if(!is.null(save_key) && variants$key[k]==save_key && !inherits(z,"error"))
     saveRDS(list(preprocessing=prep,selected_features=fs,model=fit,parameter=variants$classifier_parameter[k],variant=variants[k,]),file.path(outdir,"selected_frozen_model.rds"))
   }
  }
  cat("selector",s,sp,"completed",sum(status$status=="completed"),"/",nrow(status),"\n");flush.console()
 }
 write.table(status,file.path(outdir,"variant_status.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 saveRDS(pr,file.path(outdir,"predictions.rds"))
 list(predictions=pr,status=status)
}
