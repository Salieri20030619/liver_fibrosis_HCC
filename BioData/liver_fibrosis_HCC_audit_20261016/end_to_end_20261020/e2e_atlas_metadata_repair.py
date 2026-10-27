from pathlib import Path
import time,subprocess,os,json,traceback
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';old=P/'core55_Seurat_Harmony';O=P/'core55_Seurat_Harmony_repaired_v3';O.mkdir(exist_ok=False)
try:
 for i in range(80):
  if (old/'validation.json').exists():
   (O/'not_needed.json').write_text(json.dumps({'status':'original completed; no duplicate rerun'}));raise SystemExit(0)
  if (old/'error.txt').exists():break
  print('WAIT existing original run, no duplicate computation',i,flush=True);time.sleep(30)
 else:raise RuntimeError('Existing run still active/uncertain; no retry launched')
 text=(old/'R.log').read_text(errors='replace');assert 'No cell overlap between new meta data and Seurat object' in text,'Different error, inspect before retry'
 r=(old/'run.R').read_text();r=r.replace(str(old),str(O));r=r.replace("s$provisional_celltype<-labels[as.character(s$seurat_clusters)]","s$provisional_celltype<-setNames(unname(labels[as.character(s$seurat_clusters)]),colnames(s))")
 r=r.replace("s<-FindNeighbors(s", "saveRDS(s,'"+str(O)+"/post_Harmony_checkpoint.rds',compress=FALSE);s<-FindNeighbors(s")
 r=r.replace("markers<-FindAllMarkers", "saveRDS(s,'"+str(O)+"/annotated_checkpoint.rds',compress=FALSE);markers<-FindAllMarkers")
 r="options(error=function(){traceback(10);quit(status=1)})\n"+r
 (O/'run.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
 with open(O/'R.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run.R')],stdout=log,stderr=subprocess.STDOUT,env=env)
 assert z.returncode==0,'R repair failed; see R.log'
 import pandas as pd
 m=pd.read_csv(O/'cell_metadata_with_annotation.tsv',sep='\t',index_col=0);assert m.index.is_unique and m.donor.nunique()==55;tab=pd.crosstab(m.donor,m.provisional_celltype);tab.div(tab.sum(axis=1),axis=0).join(m.groupby('donor').F.first()).to_csv(O/'donor_celltype_fractions.tsv',sep='\t')
 v={'status':'completed','cells':len(m),'donors':55,'clusters':int(m.seurat_clusters.nunique()),'provisional_labels':m.provisional_celltype.value_counts().to_dict(),'method':'Seurat5.6 Normalize/HVG3000/Scale/PCA30/classicHarmony1.2.4 study+donor/graph0.5/UMAP; marker annotation then donor fractions','repair':'Cluster-named vector passed to Seurat metadata assignment caused zero cell-name overlap. Repaired with exact cellID-named vector; waited for original failure, no active run interrupted.','limitations':['G202379 missing MT features: partialQC','Broad marker annotations require biological review; no predefinedDTNA','No HCC samples or supplementalconflicts mixed into55F discovery','Original Seurat4.3 versus available5.6 version difference documented']};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
