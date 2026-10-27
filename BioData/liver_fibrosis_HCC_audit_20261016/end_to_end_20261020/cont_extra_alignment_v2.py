from pathlib import Path
import gzip,json,hashlib
import pandas as pd,numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;I=P/'continuation_extra_bulk_validation_v1';O=P/'continuation_extra_bulk_alignment_v2';O.mkdir(exist_ok=False)
records=[];cur={}
with gzip.open(P/'continuation_extra_bulk_metadata_v1/GSE126848_family.soft.gz','rt') as f:
 for line in f:
  if line.startswith('^SAMPLE = '):
   if cur:records.append(cur)
   cur={'GSM':line.strip().split(' = ',1)[1]}
  elif cur and line.startswith(('!Sample_title = ','!Sample_description = ')):
   k,v=line.strip().split(' = ',1);cur[k.replace('!Sample_','')]=v
  elif cur and line.startswith('!Sample_characteristics_ch1 = '):
   v=line.strip().split(' = ',1)[1]
   if ': ' in v:k,v=v.split(': ',1);cur[k]=v
 if cur:records.append(cur)
m=pd.DataFrame(records);assert m.description.str.fullmatch(r'\d+').all()
m['expression_id']=m.description.str.zfill(4);m=m.set_index('expression_id');assert m.index.is_unique and m.GSM.is_unique
x=pd.read_csv(I/'GSE126848/log2CPM_full.tsv.gz',sep='\t',index_col=0)
assert len(m)==57 and set(m.index)==set(x.columns)
m=m.loc[x.columns];m['clinical_group']=m['disease'];m.to_csv(O/'GSE126848_sample_design_verified.tsv',sep='\t')
s=pd.read_csv(I/'GSE126848/fibrosis_ssGSEA.tsv',sep='\t',index_col=0).iloc[0];assert s.index.equals(x.columns)
tfs=sorted(pd.read_csv(P/'continuation_SCENIC_candidates_v1/F5G_candidates.tsv',sep='\t').TF.unique())
rows=[]
for group,g in m.groupby('clinical_group'):
 for gene in tfs+['PDGFC']:
  for target in ['HP_HEPATIC_FIBROSIS','PDGFC','LDHB','D2HGDH']:
   if gene==target:continue
   y=s.loc[g.index] if target=='HP_HEPATIC_FIBROSIS' else x.loc[target,g.index]
   r,p=spearmanr(x.loc[gene,g.index],y);rows.append(dict(cohort='GSE126848',stratum=group,gene=gene,target=target,samples=len(g),rho=r,p=p))
c=pd.DataFrame(rows);c['BH_q']=np.nan
for _,z in c.groupby(['cohort','stratum','target']):
 ok=z.p.notna();c.loc[z.index[ok],'BH_q']=multipletests(z.loc[ok,'p'],method='fdr_bh')[1]
old=pd.read_csv(I/'extra_TF_ligand_fibrosis_correlations.tsv',sep='\t')
assert set(old.loc[old.cohort.eq('GSE126848'),'stratum'])=={'all'}
pd.concat([old,c],ignore_index=True).to_csv(O/'extra_TF_ligand_fibrosis_correlations_current.tsv',sep='\t',index=False)
v={'status':'completed','samples':57,'GSM_unique':True,'expression_IDs_exact_match':True,'alignment_rule':'GEO Sample_description is author numeric sample ID; zero-pad to four digits to exactly match all expression columns. No inferred disease-stage mapping.','clinical_group_counts':m.clinical_group.value_counts().to_dict(),'F_labels_available':False,'reused_frozen_expression_and_ssGSEA':True,'new_stratum_tests':len(c),'previous_unmatched_status_superseded':str(I/'validation.json')}
(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
