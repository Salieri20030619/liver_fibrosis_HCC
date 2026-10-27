from pathlib import Path
import json,tarfile,gzip,io,subprocess,os,traceback,itertools
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score,roc_curve,confusion_matrix
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'F6_S4_S6_F8_gene_validation';O.mkdir(exist_ok=False)
try:
 design=pd.read_csv(P/'F8_F9_real_full_cohort_v1/train_design.tsv',sep='\t',index_col=0);vectors={}
 with tarfile.open(P/'GSE135251_RAW.tar') as tar:
  for mem in tar:
   gsm=mem.name.split('_')[0]
   if not mem.isfile() or gsm not in design.index:continue
   v=pd.read_csv(io.BytesIO(gzip.decompress(tar.extractfile(mem).read())),sep='\t',header=None,index_col=0).iloc[:,0];vectors[gsm]=v[~v.index.str.startswith('__')]
 train=pd.DataFrame(vectors).loc[:,design.index];ext=pd.read_csv(D/'GSE162694_raw_counts.csv.gz',index_col=0);ed=pd.read_csv(P/'F8_F9_real_full_cohort_v1/external_design.tsv',sep='\t',index_col=0);ext=ext[ed.index]
 h=pd.read_csv(D/'HGNC_one_to_one_mapping_v1.tsv',sep='\t').set_index('ensembl_gene_id').symbol
 cohorts={};tf=pd.read_csv(P/'F5_to_F6_original_direction_v2/F5G_candidates.tsv',sep='\t').TF.unique();ligand=json.loads((P/'F6_ligand_to_F7_GSEA/validation.json').read_text())['ligand'];rows=[]
 for name,x,m in [('GSE135251',train,design),('GSE162694',ext,ed)]:
  norm=np.log2(x.div(x.sum(axis=0),axis=1)*1e6+1);ids=x.index.str.split('.').str[0];names=pd.Series(ids,index=x.index).map(h);norm.index=names.fillna(pd.Series(x.index,index=x.index));assert norm.index.is_unique
  norm.to_csv(O/(name+'_log2CPM_full.tsv.gz'),sep='\t');cohorts[name]=(norm,m)
  for group,ids2 in [('all',m.index)]+[(str(f),g.index) for f,g in m.groupby('F')]:
   for a,b in itertools.product(tf,[ligand,'LDHB','D2HGDH']):
    if a in norm.index and b in norm.index and len(ids2)>=4:
     rho,p=spearmanr(norm.loc[a,ids2],norm.loc[b,ids2]);rows.append({'cohort':name,'F_group':group,'n':len(ids2),'gene1':a,'gene2':b,'rho':rho,'p_descriptive':p})
 pd.DataFrame(rows).to_csv(O/'TF_ligand_lactate_correlations.tsv',sep='\t',index=False)
 r=f".libPaths(c('{P}/R_library','{D}/irGSEA_dependencies_v1/R_project_library_v1','{D}/trajectory_dependencies_v1/R_project_library_v1',.libPaths()));library(GSVA);library(BiocParallel);g<-strsplit(readLines('{D}/F6E_GSE164760_fibrosis_score_v1/HP_HEPATIC_FIBROSIS.gmt'),'\\t');s<-lapply(g,function(z)z[-c(1,2)]);names(s)<-sapply(g,function(z)z[1]);for(n in c('GSE135251','GSE162694')){{x<-as.matrix(read.delim(gzfile(paste0('{O}/',n,'_log2CPM_full.tsv.gz')),row.names=1,check.names=FALSE));z<-gsva(ssgseaParam(x,s,minSize=1,maxSize=Inf,normalize=TRUE),BPPARAM=SerialParam());write.table(z,paste0('{O}/',n,'_fibrosis_ssGSEA.tsv'),sep='\\t',quote=FALSE,col.names=NA)}}"
 (O/'run.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
 with open(O/'GSVA.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run.R')],stdout=log,stderr=subprocess.STDOUT,env=env)
 assert z.returncode==0;fr=[]
 for name,(norm,m) in cohorts.items():
  score=pd.read_csv(O/(name+'_fibrosis_ssGSEA.tsv'),sep='\t',index_col=0).iloc[0];assert score.index.equals(norm.columns)
  for group,ids2 in [('all',m.index)]+[(str(f),g.index) for f,g in m.groupby('F')]:
   rho,p=spearmanr(norm.loc[ligand,ids2],score[ids2]);fr.append({'cohort':name,'F_group':group,'n':len(ids2),'ligand':ligand,'rho':rho,'p_descriptive':p})
 pd.DataFrame(fr).to_csv(O/'ligand_fibrosis_score_correlations.tsv',sep='\t',index=False)
 geneauc=[];x,m=cohorts['GSE135251'];e,em=cohorts['GSE162694']
 for t in tf:
  if t not in x.index or t not in e.index:continue
  direction=1 if roc_auc_score(m.y,x.loc[t])>=.5 else -1
  geneauc.append({'gene':t,'direction_selected_in_training_only':direction,'training_AUC':roc_auc_score(m.y,direction*x.loc[t]),'external_AUC_same_direction':roc_auc_score(em.y,direction*e.loc[t])})
 pd.DataFrame(geneauc).to_csv(O/'single_gene_AUC.tsv',sep='\t',index=False)
 # Consolidate all12 algorithm families before choosing final configuration.
 a=pd.read_csv(P/'F8_F9_real_full_cohort_v1/model_metrics_with_XGBoost.tsv',sep='\t');b=pd.read_csv(P/'F8_remaining_R_families/model_metrics.tsv',sep='\t');allm=pd.concat([a,b],ignore_index=True).sort_values(['CV_mean_fold_AUC','model'],ascending=[False,True]);allm.to_csv(O/'all12_family_model_metrics.tsv',sep='\t',index=False)
 best=allm.iloc[0].to_dict();old=json.loads((P/'F8_F9_real_full_cohort_v1/validation.json').read_text())['selected']['model'];assert best['model']==old,'Global best changed; SHAP needs new model'
 pred=pd.read_csv(P/'F8_F9_real_full_cohort_v1/all_external_predictions.tsv',sep='\t');pred=pred[pred.model.eq(old)];rng=np.random.default_rng(20261020);aucs=[]
 for i in range(1000):
  ix=rng.integers(0,len(pred),len(pred));y=pred.y.to_numpy()[ix]
  if len(np.unique(y))==2:aucs.append(roc_auc_score(y,pred.probability.to_numpy()[ix]))
 v={'status':'completed','cohorts':[206,112],'gene_AUC_tests':len(geneauc),'all_model_configurations':len(allm),'algorithm_families':12,'single_models':22,'stacked_models':3,'selected_by_training_CV':best,'external_AUC_bootstrap95CI':np.quantile(aucs,[.025,.975]).tolist(),'bootstrap_unit':'expression sample; GSE162694 independent patient key unverified','S6_scope':'LDHB/D2HGDH individual gene correlations; original full lactate-gene list undefined, not invented','ligand_caution':'APOA1 selected by original broad myeloid CellChat result may reflect hepatocyte ambient RNA/uptake; no confirmed macrophage axis','original113':'Only22+3 reconstructed configurations executed; original91 ensemble identities/tuning unspecified; no claim exact113'}
 (O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
