from pathlib import Path
import time,json,pandas as pd,numpy as np,hashlib
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;O=P/'continuation_updated_bulk_associations_v1';O.mkdir(exist_ok=False);C=P/'continuation_SCENIC_candidates_v1';L=P/'continuation_ligand_to_GSEA_v1';B=P/'F6_S4_S6_F8_gene_validation';old=P/'F8_F9_real_full_cohort_v1'
for i in range(1440):
 if (C/'validation.json').exists() and (L/'validation.json').exists():break
 if (P/'continuation_SCENIC_expanded_v1/error.txt').exists():raise RuntimeError('SCENIC upstream failed')
 time.sleep(10)
else:raise RuntimeError('Upstream pending4hours')
lv=json.loads((L/'validation.json').read_text());ligand=lv.get('ligand');tfs=sorted(pd.read_csv(C/'F5G_candidates.tsv',sep='\t').TF.unique());rows=[];aucs=[];coverage=[]
for cohort,design_name in [('GSE135251','train'),('GSE162694','external')]:
 x=pd.read_csv(B/(cohort+'_log2CPM_full.tsv.gz'),sep='\t',index_col=0);m=pd.read_csv(old/(design_name+'_design.tsv'),sep='\t',index_col=0);x=x[m.index];score=pd.read_csv(B/(cohort+'_fibrosis_ssGSEA.tsv'),sep='\t',index_col=0).iloc[0].loc[m.index]
 assert 'F' in m.columns
 coverage.append({'cohort':cohort,'samples':len(m),'matched_TFs':len(set(tfs)&set(x.index)),'missing_TFs':';'.join(sorted(set(tfs)-set(x.index))),'ligand_available':ligand in x.index if ligand else False})
 for stratum,ids in [('all',m.index)]+[(str(f),g.index) for f,g in m.groupby('F')]:
  if len(ids)<5:continue
  for a in tfs+([ligand] if ligand else []):
   if a not in x.index:continue
   for b in ['HP_HEPATIC_FIBROSIS','LDHB','D2HGDH']+([ligand] if ligand else []):
    if a==b or (b!='HP_HEPATIC_FIBROSIS' and b not in x.index):continue
    y=score.loc[ids] if b=='HP_HEPATIC_FIBROSIS' else x.loc[b,ids];rho,p=spearmanr(x.loc[a,ids],y);rows.append({'cohort':cohort,'F_stratum':stratum,'gene':a,'target':b,'n':len(ids),'rho':rho,'p':p})
 for g in tfs:
  if g in x.index:aucs.append({'cohort':cohort,'gene':g,'n':len(m),'positive_F3_F4':int(m.y.sum()),'AUC_increasing_expression':roc_auc_score(m.y,x.loc[g]),'endpoint':'F3_F4_crosssectional_not_incident_HCC'})
d=pd.DataFrame(rows);d['BH_q_descriptive']=np.nan
for _,g in d.groupby(['cohort','F_stratum','target']):
 good=g.p.notna();d.loc[g.index[good],'BH_q_descriptive']=multipletests(g.loc[good,'p'],method='fdr_bh')[1]
d.to_csv(O/'new_TF_ligand_fibrosis_lactate_correlations.tsv',sep='\t',index=False);pd.DataFrame(aucs).to_csv(O/'new_TF_single_gene_AUC.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'coverage.tsv',sep='\t',index=False)
v={'status':'completed','ligand':ligand,'TF_count':len(tfs),'coverage':coverage,'candidate_sha256':hashlib.sha256((C/'F5G_candidates.tsv').read_bytes()).hexdigest(),'scope':'Two previously used staged bulk cohorts, not four or six independent validations. LDHB and D2HGDH individual correlations only, not an undisclosed complete lactate signature. F-stratified and overall descriptive results retained; no direction flipping for single-gene AUC. Previously inspected external data not pristine confirmation.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
