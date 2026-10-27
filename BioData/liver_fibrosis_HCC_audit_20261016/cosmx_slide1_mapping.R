r <- '/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016'
o <- file.path(r,'CosMx_Slide1_mapping_evidence');stopifnot(!dir.exists(o));dir.create(o)
x <- readRDS(file.path(r,'GSM8855042_Slide1_SeuratObject.RDS.complete'))
get <- function(x,k) attr(x,k,exact=TRUE)
write <- function(x,n) write.table(x,file.path(o,n),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
m <- get(x,'meta.data'); ids<-rownames(m)
write(data.frame(rds_row_1based=seq_along(ids),rds_row_id=ids,m[,intersect(c('cell_id','cell_ID','fov','orig.ident','Run_Tissue_name','slide_ID_numeric'),names(m)),drop=FALSE],check.names=FALSE),'cell_ID_fields.tsv')
sel<-grep('clusters$|assignments$|qcFlags|patient|donor|sample|tissue|slide|orig.ident',names(m),ignore.case=TRUE,value=TRUE)
write(data.frame(rds_row_id=ids,m[,unique(c('fov',sel)),drop=FALSE],check.names=FALSE),'selected_cell_annotations.tsv')
assays<-get(x,'assays'); report<-list()
for(n in names(assays)) {
 a<-assays[[n]]
 for(k in intersect(c('counts','data','scale.data'),names(attributes(a)))) {
  z<-get(a,k);dn<-get(z,'Dimnames');if(is.null(dn))dn<-dimnames(z)
  dims<-get(z,'Dim');if(is.null(dims))dims<-dim(z)
  report[[length(report)+1]]<-data.frame(assay=n,layer=k,rows=if(length(dims))dims[1] else NA,columns=if(length(dims)>1)dims[2] else NA,column_IDs=if(length(dn)>1)length(dn[[2]]) else 0,same_order=if(length(dn)>1)identical(dn[[2]],ids) else FALSE,same_set=if(length(dn)>1)setequal(dn[[2]],ids) else FALSE)
  if(length(dn)>1 && length(dn[[2]]))write(data.frame(column_1based=seq_along(dn[[2]]),cell_id=dn[[2]]),paste0('assay_',n,'_',k,'_cell_IDs.tsv'))
 }
}
if(length(report))write(do.call(rbind,report),'assay_ID_alignment.tsv')
misc<-get(x,'misc');if(is.data.frame(misc$latest.fovs))write(misc$latest.fovs,'latest_fovs.tsv')
sink(file.path(o,'mapping_structure.txt'))
cat('All metadata fields:\n');print(names(m));cat('\nProject name:\n');print(get(x,'project.name'))
for(k in c('misc','tools','images')) {
 v<-get(x,k);cat('\nATTRIBUTE',k,'length',length(v),'\n')
 for(n in names(v)) {z<-v[[n]];cat('NAME',n,'CLASS',class(z),'DIM',dim(z),'\n');if(is.data.frame(z))print(names(z));if(is.list(z)&&!is.data.frame(z))print(names(z))}
}
sink();cat('Exported annotation and assay ID evidence\n')
