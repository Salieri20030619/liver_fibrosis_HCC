from pathlib import Path
import pandas as pd,numpy as np,json
P=Path(__file__).parent;I=P/'core55_Seurat_Harmony_repaired_v3';O=P/'core55_marker_DE_annotation_v5';O.mkdir(exist_ok=False)
m=pd.read_csv(I/'cell_metadata_with_annotation.tsv',sep='\t',index_col=0);de=pd.read_csv(I/'cluster_marker_DE.tsv',sep='\t')
sets={'Hepatocyte':['ALB','APOA1','TTR','CYP2E1','CPS1','ASGR1'],'Endothelial':['PECAM1','VWF','KDR','EMCN','CLEC4G'],'Myeloid':['LYZ','LST1','TYROBP','FCER1G','CSF1R','MARCO'],'T_NK':['CD3D','CD3E','NKG7','GNLY'],'B':['MS4A1','CD79A','CD79B'],'Mesenchymal_stromal':['COL1A1','COL1A2','DCN','LUM','PDGFRA','HGF'],'Mesenchymal_mural':['RGS5','CSPG4','MCAM','NOTCH3','MYH11'],'Cholangiocyte':['KRT19','KRT7','EPCAM','SOX9']}
rows=[];labels={}
for c in sorted(m.seurat_clusters.unique()):
 v=de[de.cluster.eq(c)];passing=[]
 for label,genes in sets.items():
  q=v[v.gene.isin(genes)&v.avg_log2FC.gt(1)&v['pct.1'].ge(.2)&v['pct.1'].ge(2*v['pct.2'])]
  score=float(q.avg_log2FC.median()) if len(q)>=2 else None
  rows.append({'cluster':int(c),'candidate_label':label,'supporting_genes':';'.join(q.gene),'support_count':len(q),'median_log2FC':score})
  if score is not None:passing.append((score,label))
 passing.sort(reverse=True)
 labels[c]=passing[0][1] if passing and (len(passing)==1 or passing[0][0]-passing[1][0]>=1) else 'Unresolved_mixed'
 # Cycling genes define a state, not a validated lineage; do not assign tiny cycling cluster to a pure celltype.
 cyc=v[v.gene.isin(['MKI67','TOP2A','ANLN','BUB1B'])&v.avg_log2FC.gt(2)&v['pct.1'].gt(.2)]
 if len(cyc)>=3:labels[c]='Cycling_unresolved'
pd.DataFrame(rows).to_csv(O/'cluster_marker_support.tsv',sep='\t',index=False);pd.Series(labels,name='provisional_celltype').rename_axis('cluster').to_csv(O/'cluster_labels.tsv',sep='\t')
m['superseded_module_argmax_label']=m.provisional_celltype;m['provisional_celltype']=m.seurat_clusters.map(labels);m.to_csv(O/'cell_metadata_with_annotation.tsv',sep='\t')
tab=pd.crosstab(m.donor,m.provisional_celltype);tab.div(tab.sum(axis=1),axis=0).join(m.groupby('donor').F.first()).to_csv(O/'donor_celltype_fractions.tsv',sep='\t')
v={'status':'completed_provisional_annotation','cells':len(m),'donors':m.donor.nunique(),'clusters':m.seurat_clusters.nunique(),'labels':m.provisional_celltype.value_counts().to_dict(),'cluster_labels':{str(k):v for k,v in labels.items()},'rule':'>=2 canonical markers with log2FC>1, pct>=0.2 and >=2x other cells. Highest median marker log2FC with margin>=1 if competing labels; otherwise unresolved. Cycling>=3 markers kept lineage-unresolved. Thresholds reconstructed, not original paper parameters.','supersedes':'Only provisional labels from core55_Seurat_Harmony_repaired_v3; embeddings/graph unchanged. New annotations are sidecar TSV, original atlas.rds retains old labels.','limitations':'Marker-supported enrichment, not pure lineage, no doublet/ambient removal claimed, no predefined DTNA. DE p-values can numerically underflow, annotation uses effect and detection not exact p=0.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
