args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
src <- file.path(root, 'MASLD_snRNA_seq_seurat_v4.rds')
out <- file.path(root, 'metadata_audit')
dir.create(out, showWarnings=FALSE)
sink(file.path(out,'read_metadata.log'), split=TRUE)
cat('Start UTC:',format(Sys.time(),tz='UTC',usetz=TRUE),'\n')
cat('R version:',R.version.string,'\n')
obj <- readRDS(src)
cat('Object classes:',paste(class(obj),collapse=','),'\n')
cat('Object attribute names:',paste(names(attributes(obj)),collapse=','),'\n')
md <- attr(obj,'meta.data',exact=TRUE)
stopifnot(is.data.frame(md))
cat('Metadata dimensions:',nrow(md),ncol(md),'\n')
cat('Unique metadata rownames:',length(unique(rownames(md))),'\n')
catalog <- do.call(rbind,lapply(names(md),function(n){x<-md[[n]];data.frame(column=n,class=paste(class(x),collapse='/'),missing=sum(is.na(x)),unique_nonmissing=length(unique(x[!is.na(x)])))}))
write.csv(catalog,file.path(out,'metadata_column_catalog.csv'),row.names=FALSE)
write.csv(head(md,100),file.path(out,'metadata_first100.csv'),row.names=TRUE)
saveRDS(md,file.path(out,'metadata_only.rds'))
small <- list()
for(n in names(md)){
 u<-unique(md[[n]])
 if(length(u)<=100)small[[n]]<-data.frame(column=n,value=as.character(u),stringsAsFactors=FALSE)
}
if(length(small))write.csv(do.call(rbind,small),file.path(out,'low_cardinality_values.csv'),row.names=FALSE)
clinical <- grep('orig.ident|donor|sample|patient|FL_ID|disease|fibros|stage|SAF|NAS|biopsy|sex|gender|age|diagnos|histolog',names(md),ignore.case=TRUE,value=TRUE)
cat('Candidate clinical columns:',paste(clinical,collapse=','),'\n')
if(length(clinical))write.csv(unique(md[,clinical,drop=FALSE]),file.path(out,'clinical_distinct_rows.csv'),row.names=FALSE)
if('orig.ident'%in%names(md))write.csv(as.data.frame(table(md$orig.ident,useNA='ifany')),file.path(out,'donor_cell_counts.csv'),row.names=FALSE)
## Inspect non-expression metadata containers without executing stored commands.
tree <- list(); frame_index <- 0L; visited <- 0L
walk <- function(x,path,depth=0L){
 visited <<- visited+1L
 tree[[length(tree)+1L]] <<- data.frame(path=path,class=paste(class(x),collapse='/'),length=length(x),rows=if(is.data.frame(x))nrow(x) else NA_integer_,columns=if(is.data.frame(x))paste(names(x),collapse=';') else '',stringsAsFactors=FALSE)
 if(is.data.frame(x)){
  frame_index <<- frame_index+1L
  write.csv(head(x,10000),file.path(out,paste0('aux_metadata_frame_',frame_index,'.csv')),row.names=TRUE)
  return(invisible(NULL))
 }
 if(is.list(x)&&depth<6L&&visited<5000L){
  nn<-names(x);if(is.null(nn))nn<-as.character(seq_along(x))
  for(j in seq_along(x))walk(x[[j]],paste0(path,'/',nn[j]),depth+1L)
 }
}
for(slotname in c('misc','tools')){
 x<-attr(obj,slotname,exact=TRUE)
 if(!is.null(x))walk(x,slotname)
}
if(length(tree))write.csv(do.call(rbind,tree),file.path(out,'aux_metadata_container_catalog.csv'),row.names=FALSE)
writeLines(c(paste('metadata_rows',nrow(md),sep='\t'),paste('metadata_columns',ncol(md),sep='\t'),paste('unique_cell_ids',length(unique(rownames(md))),sep='\t'),paste('source_file',src,sep='\t')),file.path(out,'metadata_summary.tsv'))
cat('Completed UTC:',format(Sys.time(),tz='UTC',usetz=TRUE),'\n')
sink()
