from pathlib import Path
import subprocess,os,json,traceback
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';I=P/'continuation_newnetwork_models_v2'/'features';O=P/'continuation_newnetwork_models_v2'/'R_base';O.mkdir(exist_ok=False)
r=r'''.libPaths(c('LIB','OLD/irGSEA_dependencies_v1/R_project_library_v1','OLD/trajectory_dependencies_v1/R_project_library_v1',.libPaths()))
library(glmnet);library(mboost);library(plsRglm)
x<-as.matrix(read.delim('INPUT/train_features.tsv',row.names=1,check.names=FALSE));e<-as.matrix(read.delim('INPUT/external_features.tsv',row.names=1,check.names=FALSE));y<-read.delim('INPUT/train_design.tsv',row.names=1)$y;ey<-read.delim('INPUT/external_design.tsv',row.names=1)$y;fold<-read.delim('INPUT/outer_cv_folds.tsv')$fold
auc<-function(y,p){r<-rank(p);n1<-sum(y==1);n0<-sum(y==0);(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0)}
prepare<-function(a,b,c,y){mu<-colMeans(a);sdv<-apply(a,2,sd);sdv[sdv==0]<-1;a<-scale(a,mu,sdv);b<-scale(b,mu,sdv);c<-scale(c,mu,sdv);f<-colnames(a);while(length(f)>max(3,floor(ncol(a)/2))){fit<-glmnet(a[,f,drop=FALSE],y,family='binomial',alpha=0,lambda=1/nrow(a),standardize=FALSE);coefv<-abs(as.numeric(coef(fit))[-1]);f<-f[-which.min(coefv)]};list(a=as.data.frame(a[,f,drop=FALSE]),b=as.data.frame(b[,f,drop=FALSE]),c=as.data.frame(c[,f,drop=FALSE]),features=f,mean=mu,sd=sdv)}
configs<-c('Stepglm_backward','Stepglm_both','glmBoost_100','glmBoost_200','plsRglm_2','plsRglm_3');metrics<-list();external<-list();oofall<-list()
fitpredict<-function(name,a,y,b){
 if(grepl('Stepglm',name)){d<-a;d$y<-y;fit<-step(glm(y~.,data=d,family=binomial()),direction=if(grepl('backward',name))'backward' else 'both',trace=0);pred<-predict(fit,newdata=b,type='response')}
 if(grepl('glmBoost',name)){d<-a;d$y<-factor(y,levels=c(0,1));fit<-glmboost(y~.,data=d,family=Binomial(),control=boost_control(mstop=as.integer(sub('glmBoost_','',name))));pred<-predict(fit,newdata=b,type='response')}
 if(grepl('plsRglm',name)){fit<-plsRglm(object=y,dataX=a,nt=as.integer(sub('plsRglm_','',name)),modele='pls-glm-family',family='binomial',verbose=FALSE);pred<-predict(fit,newdata=b,type='response',verbose=FALSE)}
 list(model=fit,p=as.numeric(pred))
}
for(name in configs){tryCatch({
 out<-rep(NA_real_,length(y));for(k in 0:9){train<-fold!=k;test<-!train;z<-prepare(x[train,,drop=FALSE],x[test,,drop=FALSE],e,y[train]);out[test]<-fitpredict(name,z$a,y[train],z$b)$p}
 z<-prepare(x,e,e,y);fit<-fitpredict(name,z$a,y,z$b);stopifnot(all(is.finite(out)),all(is.finite(fit$p)),length(fit$p)==112)
 fa<-sapply(0:9,function(k)auc(y[fold==k],out[fold==k]));metrics[[name]]<-data.frame(model=name,CV_mean_fold_AUC=mean(fa),CV_fold_SD=sd(fa),CV_pooled_OOF_AUC=auc(y,out),external_AUC=auc(ey,fit$p),external_Brier=mean((ey-fit$p)^2));external[[name]]<-data.frame(model=name,sample=rownames(e),y=ey,probability=fit$p);oofall[[name]]<-data.frame(model=name,sample=rownames(x),fold=fold,y=y,probability=out);saveRDS(list(fit=fit$model,preprocessing=z[c('features','mean','sd')]),paste0('OUTPUT/',name,'.rds'));print(metrics[[name]])
 },error=function(err)writeLines(conditionMessage(err),paste0('OUTPUT/',name,'_error.txt')))}
write.table(do.call(rbind,metrics),'OUTPUT/model_metrics.tsv',sep='\t',row.names=FALSE,quote=FALSE);write.table(do.call(rbind,external),'OUTPUT/external_predictions.tsv',sep='\t',row.names=FALSE,quote=FALSE);write.table(do.call(rbind,oofall),'OUTPUT/training_OOF_predictions.tsv',sep='\t',row.names=FALSE,quote=FALSE);writeLines(capture.output(sessionInfo()),'OUTPUT/sessionInfo.txt')
'''.replace('LIB',str(P/'R_library')).replace('OLD',str(D)).replace('INPUT',str(I)).replace('OUTPUT',str(O))
(O/'run.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
with open(O/'R.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run.R')],stdout=log,stderr=subprocess.STDOUT,env=env)
import pandas as pd
m=pd.read_csv(O/'model_metrics.tsv',sep='\t');v={'status':'completed' if len(m)==6 else 'partial','returncode':z.returncode,'configurations':len(m),'models':m.model.tolist(),'errors':{p.name:p.read_text() for p in O.glob('*_error.txt')},'adaptation':'Same source features and outer10folds; train-fold ridge-logistic RFE, scaling train-only; native stats step, mboost glmboost and plsRglm methods; parameters reconstructed not original undisclosed grids.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
