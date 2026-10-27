from pathlib import Path
import subprocess,os,time,json,traceback
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';I=P/'mouse38_full_QC';O=P/'mouse38_Seurat_Harmony';O.mkdir(exist_ok=False)
try:
 for attempt in range(60):
  if (I/'validation.json').exists() and (P/'R_harmony_classic/harmony/DESCRIPTION').exists():break
  print('Wait',attempt,'Qc', (I/'validation.json').exists(),'harmony',(P/'R_harmony_classic/harmony/DESCRIPTION').exists(),flush=True);time.sleep(30)
 else:raise RuntimeError('Upstream QC/Harmony not ready within30min; no atlas attempted')
 r=r'''.libPaths(c('CLASSIC','LIB','OLD/irGSEA_dependencies_v1/R_project_library_v1','OLD/trajectory_dependencies_v1/R_project_library_v1',.libPaths()))
 library(Seurat);library(harmony);library(Matrix)
 set.seed(20261020)
 x<-readMM('INPUT/counts.mtx');genes<-readLines('INPUT/genes.tsv');m<-read.delim('INPUT/cell_metadata.tsv',row.names=1,check.names=FALSE)
 rownames(x)<-genes;colnames(x)<-rownames(m);stopifnot(nrow(x)==length(genes),ncol(x)==nrow(m))
 s<-CreateSeuratObject(x,meta.data=m);s<-NormalizeData(s,verbose=FALSE);s<-FindVariableFeatures(s,nfeatures=3000,selection.method='vst',verbose=FALSE);s<-ScaleData(s,features=VariableFeatures(s),verbose=FALSE);s<-RunPCA(s,npcs=30,verbose=FALSE)
 s<-RunHarmony(s,group.by.vars='FileName',reduction.use='pca',dims.use=1:30,verbose=TRUE)
 saveRDS(s,'OUTPUT/post_Harmony_checkpoint.rds',compress=FALSE);s<-FindNeighbors(s,reduction='harmony',dims=1:30,verbose=FALSE);s<-FindClusters(s,resolution=0.5,verbose=FALSE);s<-RunUMAP(s,reduction='harmony',dims=1:30,seed.use=20261020,verbose=FALSE)
 sets<-list(Hepatocyte=c('Alb','Apoa1','Ttr','Cyp2e1'),Endothelial=c('Pecam1','Vwf','Kdr','Emcn'),Myeloid=c('Lyz2','Lst1','Tyrobp','Fcer1g'),T_NK=c('Cd3d','Cd3e','Nkg7','Gnly'),B=c('Ms4a1','Cd79a','Cd79b'),Mesenchymal_stromal=c('Col1a1','Col1a2','Dcn','Lum','Pdgfra'),Mesenchymal_mural=c('Rgs5','Cspg4','Mcam','Notch3'),Cholangiocyte=c('Krt19','Krt7','Epcam','Sox9'))
 sets<-lapply(sets,function(z)intersect(z,rownames(s)));stopifnot(all(lengths(sets)>=2));s<-AddModuleScore(s,features=sets,name='markerScore',seed=20261020)
 scores<-s@meta.data[,paste0('markerScore',seq_along(sets))];colnames(scores)<-names(sets);agg<-aggregate(scores,list(cluster=s$seurat_clusters),mean);labels<-apply(agg[,-1],1,function(v){o<-order(v,decreasing=TRUE);if(v[o[1]]<=0 || v[o[1]]-v[o[2]]<0.1)'Unresolved_mixed' else names(sets)[o[1]]});names(labels)<-agg$cluster;s$provisional_celltype<-setNames(unname(labels[as.character(s$seurat_clusters)]),colnames(s))
 write.table(agg,'OUTPUT/cluster_marker_scores.tsv',sep='\t',quote=FALSE,row.names=FALSE);write.table(s@meta.data,'OUTPUT/cell_metadata_with_annotation.tsv',sep='\t',quote=FALSE,col.names=NA)
 markers<-FindAllMarkers(s,only.pos=TRUE,test.use='wilcox',min.pct=.25,logfc.threshold=.25,verbose=FALSE);write.table(markers,'OUTPUT/cluster_marker_DE.tsv',sep='\t',quote=FALSE,row.names=FALSE)
 write.table(Embeddings(s,'harmony'),'OUTPUT/harmony_coordinates.tsv',sep='\t',quote=FALSE,col.names=NA);write.table(Embeddings(s,'umap'),'OUTPUT/umap_coordinates.tsv',sep='\t',quote=FALSE,col.names=NA)
 saveRDS(s,'OUTPUT/atlas.rds');writeLines(capture.output(sessionInfo()),'OUTPUT/sessionInfo.txt');writeLines('completed','OUTPUT/R_COMPLETE.txt')
 '''.replace('CLASSIC',str(P/'R_harmony_classic')).replace('LIB',str(P/'R_library')).replace('OLD',str(D)).replace('INPUT',str(I)).replace('OUTPUT',str(O))
 (O/'run.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
 with open(O/'R.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run.R')],stdout=log,stderr=subprocess.STDOUT,env=env)
 assert z.returncode==0,'R execution failed; see R.log'
 import pandas as pd,numpy as np
 m=pd.read_csv(O/'cell_metadata_with_annotation.tsv',sep='\t',index_col=0);assert m.index.is_unique and m.FileName.nunique()==38
 # Downstream donor-level abundance reads the persisted annotation.
 tab=pd.crosstab(m.donor,m.provisional_celltype);frac=tab.div(tab.sum(axis=1),axis=0);frac.join(m.groupby('donor').F.first()).to_csv(O/'donor_celltype_fractions.tsv',sep='\t')
 (O/'validation.json').write_text(json.dumps({'status':'completed','cells':len(m),'source_units':int(m.donor.nunique()),'libraries':38,'clusters':int(m.seurat_clusters.nunique()),'provisional_labels':m.provisional_celltype.value_counts().to_dict(),'method':'Seurat Normalize/HVG3000/Scale/PCA30/Harmony donor/graph resolution0.5/UMAP; parameters beyond original specification reconstructed','version_deviation':'Seurat5.6/R4.5 versus original4.3; classicHarmony1.2.4 used explicitly, notHarmony2. Missing MT features prevent full MTQC; numeric pipeline can proceed with flagged partialQC','annotation_limit':'Marker-score coarse labels provisional; no DTNA or axis assumed. Not a validated pure-HSC reference; numeric clusters require review.','integration_limit':'Full38 mouse libraries only; capture/timepoint confounding can remove biological differences. Sourceunits are not asserted independentanimals. Author celltype labels retained separately.'},indent=2));print('Complete',len(m),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
