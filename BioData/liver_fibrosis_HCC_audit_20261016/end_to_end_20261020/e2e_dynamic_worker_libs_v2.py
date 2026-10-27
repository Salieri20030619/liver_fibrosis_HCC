from pathlib import Path
import os,subprocess,json,time
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'F3_F5_dynamic_genes_workerlibs_v2';O.mkdir(exist_ok=False)
libs=[P/'R_monocle_legacy',P/'R_library',D/'trajectory_dependencies_v1/R_project_library_v1',D/'irGSEA_dependencies_v1/R_project_library_v1']
env=os.environ.copy();env['R_LIBS']=':'.join(map(str,libs));n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
records=[]
for branch in ['F5','F3']:
 src=P/'F3_F5_dynamic_genes_resumed';out=O/branch
 for i in range(100):
  v=json.loads((src/'validation.json').read_text()) if (src/'validation.json').exists() else {}
  prior=next((x for x in v.get('branches',[]) if x['branch']==branch),None)
  if prior:break
  print('WAIT preceding branch completion',branch,i,flush=True);time.sleep(30)
 if not prior:raise RuntimeError('Prior branch not finished; no duplicate started')
 out.mkdir()
 if prior['status']=='completed':records.append(prior);continue
 checkpoint=src/branch/'dispersion_checkpoint.rds'
 if not checkpoint.exists():records.append({'branch':branch,'status':'blocked','reason':'No dispersion checkpoint'});continue
 r=f"library(monocle);stopifnot(as.character(packageVersion('dplyr'))=='1.1.4');cds<-readRDS('{checkpoint}');genes<-rownames(fData(cds))[fData(cds)$use_for_ordering];writeLines(as.character(length(genes)),'{out}/tested_gene_count.txt');de<-differentialGeneTest(cds[genes,],fullModelFormulaStr='~sm.ns(Pseudotime,df=3)',reducedModelFormulaStr='~1',cores=4);write.table(de,'{out}/pseudotime_dynamic_genes.tsv',sep='\\t',quote=FALSE,col.names=NA);writeLines(capture.output(sessionInfo()),'{out}/sessionInfo.txt');writeLines('completed','{out}/COMPLETE.txt')"
 (out/'run.R').write_text(r)
 with open(out/'R.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(out/'run.R')],env=env,stdout=log,stderr=subprocess.STDOUT)
 records.append({'branch':branch,'status':'completed' if (out/'COMPLETE.txt').exists() else 'failed','returncode':z.returncode,'checkpoint':str(checkpoint)})
 (O/'validation.json').write_text(json.dumps({'status':'running' if branch=='F5' else 'finished','branches':records,'repair':'R_LIBS inherited by PSOCK workers, same Monocle tests, resume from dispersion checkpoint'},indent=2));print(json.dumps(records[-1]),flush=True)
