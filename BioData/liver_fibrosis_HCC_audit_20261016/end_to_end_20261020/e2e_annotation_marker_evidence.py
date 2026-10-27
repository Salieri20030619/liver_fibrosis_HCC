from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'core55_annotation_marker_review_v4';O.mkdir(exist_ok=False)
m=pd.read_csv(P/'core55_Seurat_Harmony_repaired_v3/cell_metadata_with_annotation.tsv',sep='\t',index_col=0)
sets={'Hepatocyte':['ALB','APOA1','TTR','CYP2E1'],'Endothelial':['PECAM1','VWF','KDR','EMCN'],'Myeloid':['LYZ','LST1','TYROBP','FCER1G'],'T_NK':['CD3D','CD3E','NKG7','GNLY'],'B':['MS4A1','CD79A','CD79B'],'Mesenchymal_stromal':['COL1A1','COL1A2','DCN','LUM','PDGFRA'],'Mesenchymal_mural':['RGS5','CSPG4','MCAM','NOTCH3'],'Cholangiocyte':['KRT19','KRT7','EPCAM','SOX9']};genes=sum(sets.values(),[]);rows=[]
for study,file in [('GSE202379',P/'GSE202379_core39_QC/core39_counts.h5ad'),('GSE244832',D/'GSE244832/qc_pass_counts.h5ad')]:
 a=ad.read_h5ad(file);a.obs_names=study+'|'+a.obs_names
 if study=='GSE244832':
  sy=a.var.gene_symbol_original.astype(str);a=a[:,~sy.duplicated(keep=False)].copy();a.var_names=a.var.gene_symbol_original.astype(str)
 ids=m.index[m.study.eq(study)];a=a[ids];total=np.asarray(a.X.sum(axis=1)).ravel();g=[x for x in genes if x in a.var_names];raw=a[:,g].X.toarray();v=np.log1p(raw/np.maximum(total[:,None],1)*1e4);clusters=m.loc[ids,'seurat_clusters'].to_numpy()
 for c in np.unique(clusters):
  keep=clusters==c
  for j,gene in enumerate(g):rows.append({'study':study,'cluster':int(c),'gene':gene,'cells':int(keep.sum()),'mean_log1p_CP10k':float(v[keep,j].mean()),'pct_detected':float((raw[keep,j]>0).mean())})
z=pd.DataFrame(rows);z.to_csv(O/'canonical_marker_expression_by_study_cluster.tsv',sep='\t',index=False)
out=[]
for (c,g),v in z.groupby(['cluster','gene']):out.append({'cluster':c,'gene':g,'mean_log1p_CP10k':np.average(v.mean_log1p_CP10k,weights=v.cells),'pct_detected':np.average(v.pct_detected,weights=v.cells)})
w=pd.DataFrame(out);w.to_csv(O/'canonical_marker_expression_pooled.tsv',sep='\t',index=False);(O/'validation.json').write_text(json.dumps({'status':'completed','cells':len(m),'clusters':m.seurat_clusters.nunique(),'markers':len(genes),'reason':'Different gene-set AddModuleScore values are not directly interchangeable. Strong hepatocyte signal masks otherwise specific nonparenchymal markers. Preserve original labels; inspect marker evidence before downstream.'},indent=2));print(w[w.cluster.isin([5,7,10,12,18])].to_json(orient='records'))
