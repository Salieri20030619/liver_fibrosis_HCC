from pathlib import Path
import os,subprocess,json,sys,traceback
P=Path(__file__).parent;mode=sys.argv[1];O=P/('resume_'+mode+'_v7');O.mkdir(exist_ok=False)
def run(args,env):
 with open(O/'execution.log','a',buffering=1) as f:return subprocess.run(args,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
try:
 env=os.environ.copy();env['RETICULATE_PYTHON']='/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python';native=P/'dependency_irgsea_repair_v3/native_dependencies';env['LD_PRELOAD']=str(native/'lib/libstdc++.so.6')+':'+str(native/'lib/libdeflate.so.0');env['LD_LIBRARY_PATH']=str(native/'lib');lib=P/'trajectory_dependencies_v1/R_project_library_v1';ilib=P/'irGSEA_dependencies_v1/R_project_library_v1'
 if mode=='cytotrace':
  prefix=O/'netcdf_native';env['CONDA_PKGS_DIRS']=str(O/'conda_cache');rc=run(['/mnt/BioDisk/wangmingyuan/miniconda3/bin/conda','create','--yes','--prefix',str(prefix),'--override-channels','-c','conda-forge','libnetcdf','pkg-config'],env);assert rc==0
  env['PATH']=str(prefix/'bin')+':'+env['PATH'];env['LD_LIBRARY_PATH']=str(prefix/'lib')+':'+env['LD_LIBRARY_PATH'];env['PKG_CONFIG_PATH']=str(prefix/'lib/pkgconfig')
  r=f".libPaths(c('{lib}','{ilib}',.Library));options(timeout=600,Ncpus=4,repos=BiocManager::repositories());BiocManager::install('sva',lib='{lib}',update=FALSE,ask=FALSE);install.packages('ncdf4',lib='{lib}',configure.args='--with-nc-config={prefix/'bin/nc-config'}');install.packages('HiClimR',lib='{lib}');remotes::install_local('{P/'dependency_cytotrace_code_v5/CytoTRACE_code_only.tar.gz'}',lib='{lib}',upgrade='never',dependencies=FALSE);stopifnot(requireNamespace('CytoTRACE',quietly=TRUE));writeLines(as.character(packageVersion('CytoTRACE')),'{O/'installed_version.txt'}')"
 else:
  r=f".libPaths(c('{ilib}','{lib}',.Library));options(timeout=600,repos=c(CRAN='https://cloud.r-project.org'));install.packages('doMC',lib='{ilib}');stopifnot(requireNamespace('doMC',quietly=TRUE));writeLines(as.character(packageVersion('irGSEA')),'{O/'installed_version.txt'}')"
 q=O/'install.R';q.write_text(r);rc=run(['/usr/bin/Rscript','--vanilla',str(q)],env);assert rc==0,(mode,'install',rc)
 s=(P/'stage1_original_methods_v1.py').read_text().replace('F3_CytoTRACE_official_v1','F3_CytoTRACE_official_v2').replace('F4_irGSEA_official_v1','F4_irGSEA_official_v2')
 s=s.replace("dep=P/('dependency_cytotrace_code_v5' if mode=='cytotrace' else 'dependency_irgsea_repair_v5')","dep=P/('resume_cytotrace_v7' if mode=='cytotrace' else 'resume_irgsea_v7')")
 old="obj<-irGSEA::irGSEA.score(object=x,min.cells=0,min.feature=0,seeds=0,ncores=4,custom=TRUE,geneset=sets,msigdb=FALSE,method=meth,kcdf='Gaussian',minGSSize=1,maxGSSize=Inf);"
 new="old<-'{P/'F4_irGSEA_official_v1'}/'+g"
 replacement="old<-file.path('{P/'F4_irGSEA_official_v1'}','{g}','irGSEA_cell_scores.rds');if(file.exists(old)){{obj<-readRDS(old);missing<-setdiff(meth,Assays(obj));if(length(missing)>0)obj<-irGSEA::irGSEA.score(object=obj,assay='RNA',slot='data',add=TRUE,seeds=0,ncores=4,custom=TRUE,geneset=sets,msigdb=FALSE,method=missing,kcdf='Gaussian',minGSSize=1,maxGSSize=Inf)}}else{{obj<-irGSEA::irGSEA.score(object=x,min.cells=0,min.feature=0,seeds=0,ncores=4,custom=TRUE,geneset=sets,msigdb=FALSE,method=meth,kcdf='Gaussian',minGSSize=1,maxGSSize=Inf)}};"
 assert old in s;s=s.replace(old,replacement)
 q=O/'analysis_resume.py';q.write_text(s.replace('P=Path(__file__).parent;','P=Path(__file__).parent.parent;'));compile(q.read_text(),str(q),'exec');rc=run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python',str(q),mode],env);(O/'result.json').write_text(json.dumps({'analysis_returncode':rc}))
except Exception:(O/'error.txt').write_text(traceback.format_exc());raise
