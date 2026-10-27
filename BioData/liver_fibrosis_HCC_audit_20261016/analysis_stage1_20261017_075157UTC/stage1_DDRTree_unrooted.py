import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,subprocess,traceback,numpy as np,pandas as pd,anndata as ad,scipy.io
P=Path(__file__).parent;O=P/'F3_DDRTree_unrooted_v1';O.mkdir(exist_ok=False);results=[]
for g in ['GSE244832','GSE136103']:
 try:
  d=O/g;d.mkdir();source=P/'F2_Roe_F3_F4_subclusters_v1'/(g+'_Mesenchymal');a=ad.read_h5ad(source/'subcluster_checkpoint.h5ad',backed='r');obs=a.obs.copy();genes=a.var_names.tolist();a.file.close();src=ad.read_h5ad(P/g/'qc_pass_counts.h5ad',backed='r');counts=src[obs.index,:].to_memory();src.file.close();scipy.io.mmwrite(str(d/'counts.mtx'),counts.X.T);obs.to_csv(d/'cell_metadata.tsv',sep='\t');pd.DataFrame({'gene_id':counts.var_names,'gene_short_name':counts.var.gene_symbol_original.to_numpy()}).to_csv(d/'gene_metadata.tsv',sep='\t',index=False);(d/'ordering_genes.txt').write_text('\n'.join(genes));del counts
  r="""args<-commandArgs(trailingOnly=TRUE);d<-args[1];lib<-args[2];.libPaths(c(lib,.Library));library(monocle);library(Matrix);set.seed(0)
x<-readMM(file.path(d,'counts.mtx'));cell<-read.delim(file.path(d,'cell_metadata.tsv'),row.names=1,check.names=FALSE);gene<-read.delim(file.path(d,'gene_metadata.tsv'),check.names=FALSE);rownames(gene)<-gene$gene_id;rownames(x)<-gene$gene_id;colnames(x)<-rownames(cell)
cds<-newCellDataSet(as(x,'dgCMatrix'),phenoData=new('AnnotatedDataFrame',data=cell),featureData=new('AnnotatedDataFrame',data=gene),expressionFamily=negbinomial.size())
cds<-estimateSizeFactors(cds);cds<-setOrderingFilter(cds,readLines(file.path(d,'ordering_genes.txt')))
cds<-reduceDimension(cds,max_components=2,reduction_method='DDRTree',norm_method='log',pseudo_expr=1,relative_expr=TRUE,verbose=TRUE)
saveRDS(cds,file.path(d,'DDRTree_unrooted.rds'));xy<-t(reducedDimS(cds));write.table(cbind(cell_id=rownames(xy),xy),file.path(d,'DDRTree_coordinates.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
png(file.path(d,'F3H_unrooted_DDRTree.png'),width=1500,height=1000,res=150);par(mar=c(5,5,4,10));cl<-as.factor(cell$subcluster_preliminary);plot(xy[,1],xy[,2],col=as.integer(cl),pch=16,cex=.3,xlab='DDRTree 1',ylab='DDRTree 2',main='Unrooted mesenchymal topology: no developmental direction');legend('topright',legend=levels(cl),col=seq_along(levels(cl)),pch=16,cex=.6);dev.off();writeLines(capture.output(sessionInfo()),file.path(d,'sessionInfo.txt'))
"""; (d/'run.R').write_text(r)
  with open(d/'R.log','w') as log:ret=subprocess.run(['/usr/bin/Rscript','--vanilla',str(d/'run.R'),str(d),str(P/'trajectory_dependencies_v1/R_project_library_v1')],stdout=log,stderr=subprocess.STDOUT)
  results.append({'cohort':g,'cells':len(obs),'returncode':ret.returncode,'complete':(d/'DDRTree_unrooted.rds').exists(),'ordering_genes':len(genes),'root_assigned':False});(O/'state.json').write_text(json.dumps(results));print(json.dumps(results[-1]),flush=True)
 except Exception:results.append({'cohort':g,'error':traceback.format_exc()});(O/'state.json').write_text(json.dumps(results));print(traceback.format_exc(),flush=True)
(O/'scope.json').write_text(json.dumps({'panel':'F3H partial: original Monocle2 DDRTree topology, no ordered pseudotime','input':'all original broad mesenchymal QC cells including D9/D10 only for unsupervised topology; no disease inference or pooling studies','ordering_genes':'reuse existing 3000 HVGs, implementation choice because original ordering-feature specification unclear','not_completed':'No orderCells root, CytoTRACE, directed branch assignment or dynamic-gene claims. Broad mesenchymal identity not pure HSC.'},indent=2))
