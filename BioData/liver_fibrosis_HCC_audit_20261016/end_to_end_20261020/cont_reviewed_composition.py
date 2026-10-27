from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import spearmanr,mannwhitneyu
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;O=P/'continuation_reviewed_composition_v1';O.mkdir(exist_ok=False)
m=pd.read_csv(P/'continuation_annotation_refined_v1/cell_metadata_reviewed.tsv.gz',sep='\t',index_col=0)
meta=m.groupby('donor')[['study','F']].first();assert len(meta)==55
alltotal=m.groupby('donor').size();tables=[];tests=[]
for scope,sub,field in [('all_broad',m,'reviewed_broad_lineage'),('myeloid_states',m[m.reviewed_broad_lineage.eq('Myeloid')],'refined_cluster'),('stromal_states',m[m.reviewed_broad_lineage.eq('Stromal')],'refined_cluster'),('mural_states',m[m.reviewed_broad_lineage.eq('Mural')],'refined_cluster')]:
 counts=pd.crosstab(sub.donor,sub[field]).reindex(meta.index,fill_value=0);within=counts.sum(axis=1)
 for label in counts:
  t=meta.copy();t['scope']=scope;t['label']=label;t['cells']=counts[label];t['all_QC_cells']=alltotal;t['scope_cells']=within;t['fraction_all_cells']=t.cells/t.all_QC_cells;t['fraction_within_scope']=t.cells/t.scope_cells.replace(0,np.nan);tables.append(t.reset_index())
  for study,ids in [('pooled_descriptive',meta.index)]+[(s,g.index) for s,g in meta.groupby('study')]:
   for denominator in ['fraction_all_cells','fraction_within_scope']:
    u=t.loc[ids].dropna(subset=[denominator]);rho,p=spearmanr(u.F,u[denominator]) if len(u)>3 else (np.nan,np.nan)
    tests.append(dict(scope=scope,label=label,study=study,denominator=denominator,test='F_spearman',records=len(u),rho=rho,p=p,n_high=np.nan,n_low=np.nan,median_difference=np.nan))
    for contrast,high in [('F3_F4_vs_F0_F1_F2',u.F.ge(3)),('F4_vs_F0_F1_F2_F3',u.F.eq(4))]:
     a=u.loc[high,denominator];b=u.loc[~high,denominator];valid=min(len(a),len(b))>=3
     tests.append(dict(scope=scope,label=label,study=study,denominator=denominator,test=contrast,records=len(u),rho=np.nan,p=mannwhitneyu(a,b,alternative='two-sided').pvalue if valid else np.nan,n_high=len(a),n_low=len(b),median_difference=a.median()-b.median() if valid else np.nan))
pd.concat(tables,ignore_index=True).to_csv(O/'source_record_composition.tsv',sep='\t',index=False)
d=pd.DataFrame(tests);d['BH_q_descriptive']=np.nan
for _,g in d.groupby(['scope','study','denominator','test']):
 good=g.p.notna();d.loc[g.index[good],'BH_q_descriptive']=multipletests(g.loc[good,'p'],method='fdr_bh')[1]
d.to_csv(O/'within_study_F_composition_associations.tsv',sep='\t',index=False)
v={'status':'completed','source_records':len(meta),'cells':len(m),'tests':len(d),'mixed_cells_retained_in_allcell_denominator':int(m.reviewed_broad_lineage.str.contains('mixed').sum()),'scope':'Final reviewed annotation; record-level fractions with zero counts included; no cell-level pseudoreplication. Within-scope fractions undefined for records with zero cells in scope. BH within scope/study/denominator/contrast.','limits':'Source records not proven globally independent people. Observational cross-sectional composition and capture biases; source and F1/F2 confounded. HCC excluded from core rather than merged with source. No biological progression or causal claim.'}
(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
