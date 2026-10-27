lib<-'/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/analysis_stage1_20261017_075157UTC/irGSEA_dependencies_v1/R_project_library_v1'
.libPaths(c(lib,'/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/analysis_stage1_20261017_075157UTC/trajectory_dependencies_v1/R_project_library_v1',.Library))
options(timeout=600,Ncpus=4,repos=BiocManager::repositories())
install.packages('https://cran.r-project.org/src/contrib/Archive/gghalves/gghalves_0.1.4.tar.gz',repos=NULL,type='source',lib=lib)
BiocManager::install(c('magick','SpatialExperiment','GSVA'),lib=lib,ask=FALSE,update=FALSE,Ncpus=4)
remotes::install_github('chuiqin/irGSEA',lib=lib,upgrade='never',dependencies=NA)
stopifnot(requireNamespace('irGSEA',quietly=TRUE))
writeLines(as.character(packageVersion('irGSEA')),'/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/analysis_stage1_20261017_075157UTC/dependency_irgsea_repair_v5/installed_version.txt')
