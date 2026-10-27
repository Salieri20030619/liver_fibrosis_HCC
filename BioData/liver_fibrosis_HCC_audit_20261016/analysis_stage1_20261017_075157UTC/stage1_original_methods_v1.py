from pathlib import Path
import os,sys,time,subprocess,json,traceback
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
P=Path(__file__).parent;mode=sys.argv[1];O=P/('F3_CytoTRACE_official_v1' if mode=='cytotrace' else 'F4_irGSEA_official_v1');O.mkdir(exist_ok=False)
dep=P/('dependency_cytotrace_code_v5' if mode=='cytotrace' else 'dependency_irgsea_repair_v5');pid=318035 if mode=='cytotrace' else 317062
(O/'state.json').write_text(json.dumps({'status':'waiting_for_verified_install','dependency':str(dep),'installer_pid':pid}))
try:
 for _ in range(240):
  if (dep/'installed_version.txt').exists():break
  if not Path('/proc/'+str(pid)).exists():raise RuntimeError('Installer exited without installed_version confirmation; inspect '+str(dep/'execution.log'))
  time.sleep(15)
 else:raise TimeoutError('Installation not confirmed after 1 hour')
 env=os.environ.copy();env['RETICULATE_PYTHON']='/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python';env['RETICULATE_AUTOCONFIGURE']='FALSE'
 if mode=='irgsea':
  prefix=P/'dependency_irgsea_repair_v3/native_dependencies';env['LD_PRELOAD']=str(prefix/'lib/libstdc++.so.6')+':'+str(prefix/'lib/libdeflate.so.0');env['LD_LIBRARY_PATH']=str(prefix/'lib')+':'+env.get('LD_LIBRARY_PATH','')
 for g in ['GSE136103','GSE244832']:
  d=O/g;d.mkdir();(O/'state.json').write_text(json.dumps({'status':'running','cohort':g,'method':mode}))
  if mode=='cytotrace':
   src=P/'F3_DDRTree_unrooted_v1'/g
   r=f""".libPaths(c('{P/'trajectory_dependencies_v1/R_project_library_v1'}',.Library));library(Matrix);set.seed(0)
stopifnot(requireNamespace('CytoTRACE',quietly=TRUE));src<-'{src}';out<-'{d}'
x<-as.matrix(readMM(file.path(src,'counts.mtx')));genes<-read.delim(file.path(src,'gene_metadata.tsv'),check.names=FALSE);meta<-read.delim(file.path(src,'cell_metadata.tsv'),row.names=1,check.names=FALSE);rownames(x)<-genes$gene_id;colnames(x)<-rownames(meta);stopifnot(all(is.finite(x)),all(x>=0),all(x==round(x)))
z<-CytoTRACE::CytoTRACE(x,enableFast=TRUE,ncores=4,subsamplesize=1000);saveRDS(z,file.path(out,'CytoTRACE_original_result.rds'))
tab<-data.frame(cell_id=names(z$CytoTRACE),CytoTRACE=z$CytoTRACE,detected_genes=z$Counts,GCS=z$GCS);write.table(tab,file.path(out,'cell_scores.tsv'),sep='\t',row.names=FALSE,quote=FALSE);writeLines(z$filteredCells,file.path(out,'algorithm_filtered_cells.txt'));writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
"""
  else:
   src=P/'F4_AUCell_official_v1'/g
   r=f""".libPaths(c('{P/'irGSEA_dependencies_v1/R_project_library_v1'}','{P/'trajectory_dependencies_v1/R_project_library_v1'}',.Library));library(irGSEA);library(Seurat);library(Matrix);set.seed(0)
src<-'{src}';out<-'{d}';x<-as(readMM(file.path(src,'expression.mtx')),'CsparseMatrix');rownames(x)<-readLines(file.path(src,'genes.txt'));colnames(x)<-readLines(file.path(src,'cells.txt'));x@x<-log1p(x@x);meta<-read.delim(file.path(src,'cell_metadata.tsv'),row.names=1,check.names=FALSE);stopifnot(identical(colnames(x),rownames(meta)))
gmt<-strsplit(readLines('{P/'F4_official_sources_v1/Hallmark_2025.1_Hs.gmt'}'),'\t');sets<-lapply(gmt,function(z)unique(z[-c(1,2)]));names(sets)<-sapply(gmt,function(z)z[1]);meth<-c('AUCell','UCell','singscore','ssgsea')
obj<-irGSEA::irGSEA.score(object=x,min.cells=0,min.feature=0,seeds=0,ncores=4,custom=TRUE,geneset=sets,msigdb=FALSE,method=meth,kcdf='Gaussian',minGSSize=1,maxGSSize=Inf);saveRDS(obj,file.path(out,'irGSEA_cell_scores.rds'));stopifnot(all(meth %in% Assays(obj)))
for(sort in unique(meta$sorting)){{
 m<-meta[meta$sorting==sort,,drop=FALSE];ids<-unique(m$donor);groups<-sapply(ids,function(id){{v<-unique(m$group[m$donor==id]);stopifnot(length(v)==1);v}});dm<-data.frame(group=groups,row.names=ids);tiny<-matrix(1,nrow=2,ncol=length(ids),dimnames=list(c('placeholder1','placeholder2'),ids));don<-CreateSeuratObject(counts=tiny,meta.data=dm,min.cells=0,min.features=0)
 for(method in meth){{
  scores<-GetAssayData(obj,assay=method,layer='scale.data');stopifnot(all(rownames(m) %in% colnames(scores)));avg<-sapply(ids,function(id)rowMeans(scores[,rownames(m)[m$donor==id],drop=FALSE]));rownames(avg)<-rownames(scores);colnames(avg)<-ids;stopifnot(all(is.finite(avg)));a<-CreateAssayObject(data=avg);a<-SetAssayData(a,layer='scale.data',new.data=avg);fm<-obj[[method]][[]];a<-AddMetaData(a,metadata=fm);don[[method]]<-a;write.table(avg,file.path(out,paste0(sort,'_',method,'_donor_scores.tsv')),sep='\t',quote=FALSE,col.names=NA)
 }}
 saveRDS(don,file.path(out,paste0(sort,'_donor_object.rds')));res<-tryCatch(irGSEA::irGSEA.integrate(don,group.by='group',method=meth,p.value='p_val_adj'),error=function(e){{writeLines(conditionMessage(e),file.path(out,paste0(sort,'_integration_error.txt')));NULL}})
 if(!is.null(res)){{saveRDS(res,file.path(out,paste0(sort,'_irGSEA_integration.rds')));for(n in names(res))write.table(res[[n]],file.path(out,paste0(sort,'_',n,'_integration.tsv')),sep='\t',quote=FALSE,row.names=FALSE)}}
}}
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
"""
  q=d/'run.R';q.write_text(r)
  with open(d/'R.log','w') as log:ret=subprocess.run(['/usr/bin/Rscript','--vanilla',str(q)],stdout=log,stderr=subprocess.STDOUT,env=env)
  (d/'result.json').write_text(json.dumps({'R_returncode':ret.returncode}));assert ret.returncode==0,(g,ret.returncode)
 (O/'state.json').write_text(json.dumps({'status':'R_steps_completed','method':mode,'plots':'pending validation and rendering'}))
except Exception:(O/'error.txt').write_text(traceback.format_exc());(O/'state.json').write_text(json.dumps({'status':'blocked','method':mode}));raise
