from pathlib import Path
import subprocess,os,json,time
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'F3_F5_dynamic_genes_resumed';O.mkdir(exist_ok=False)
for i in range(60):
 if (P/'R_monocle_legacy/dplyr/DESCRIPTION').exists() and not (P/'R_monocle_legacy/00LOCK-dplyr').exists():break
 print('WAIT legacy dplyr',i,flush=True);time.sleep(30)
else:raise RuntimeError('Legacy dplyr unavailable')
env=os.environ.copy();env['RETICULATE_PYTHON']='/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python';env['RETICULATE_AUTOCONFIGURE']='FALSE';n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
records=[]
for branch in ['F5','F3']:
 src=P/(branch+'_Monocle_API_repaired_v4');out=O/branch;out.mkdir()
 for i in range(80):
  if (src/'validation.json').exists():break
  print('WAIT existing',branch,'rooted checkpoint',i,flush=True);time.sleep(30)
 if not (src/'ordered_checkpoint.rds').exists():records.append({'branch':branch,'status':'blocked','reason':'No ordered checkpoint; see upstream error'});continue
 r=f".libPaths(c('{P}/R_monocle_legacy','{P}/R_library','{D}/trajectory_dependencies_v1/R_project_library_v1','{D}/irGSEA_dependencies_v1/R_project_library_v1',.libPaths()));library(monocle);stopifnot(as.character(packageVersion('dplyr'))=='1.1.4');cds<-readRDS('{src}/ordered_checkpoint.rds');cds<-estimateDispersions(cds);saveRDS(cds,'{out}/dispersion_checkpoint.rds');genes<-rownames(fData(cds))[fData(cds)$use_for_ordering];de<-differentialGeneTest(cds[genes,],fullModelFormulaStr='~sm.ns(Pseudotime,df=3)',reducedModelFormulaStr='~1',cores=4);write.table(de,'{out}/pseudotime_dynamic_genes.tsv',sep='\\t',quote=FALSE,col.names=NA);writeLines(capture.output(sessionInfo()),'{out}/sessionInfo.txt');writeLines('completed','{out}/COMPLETE.txt')"
 (out/'run.R').write_text(r)
 with open(out/'R.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(out/'run.R')],env=env,stdout=log,stderr=subprocess.STDOUT)
 records.append({'branch':branch,'status':'completed' if (out/'COMPLETE.txt').exists() else 'failed','returncode':z.returncode,'source_ordered_checkpoint':str(src/'ordered_checkpoint.rds')});(O/'validation.json').write_text(json.dumps({'status':'running' if branch=='F5' else 'finished','branches':records,'repair':'Legacydplyr1.1.4 projectonly for Monocle group_by_ compatibility; no trajectory recomputation','inference_limit':'Exploratory cell-level pseudotime associations, not patient progression or causality'},indent=2));print(json.dumps(records[-1]),flush=True)
