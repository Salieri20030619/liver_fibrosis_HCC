from pathlib import Path
import json,time,hashlib,os,subprocess,gzip,re,numpy as np,pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';I=P/'continuation_extra_bulk_inputs_v1';O=P/'continuation_extra_bulk_validation_v1';O.mkdir(exist_ok=False);tfpath=P/'continuation_SCENIC_candidates_v1/F5G_candidates.tsv';tfs=sorted(pd.read_csv(tfpath,sep='\t').TF.unique());ligand=json.loads((P/'continuation_ligand_to_GSEA_v1/validation.json').read_text())['ligand'];rows=[];corr=[];h=pd.read_csv(D/'HGNC_one_to_one_mapping_v1.tsv',sep='\t').set_index('ensembl_gene_id').symbol
def soft_meta(acc):
 records=[];cur={}
 with gzip.open(P/'continuation_extra_bulk_metadata_v1'/(acc+'_family.soft.gz'),'rt') as f:
  for line in f:
   if line.startswith('^SAMPLE = '):
    if cur:records.append(cur)
    cur={'GSM':line.strip().split(' = ',1)[1]}
   elif cur and line.startswith('!Sample_title = '):cur['title']=line.strip().split(' = ',1)[1]
   elif cur and line.startswith('!Sample_characteristics_ch1 = '):
    v=line.strip().split(' = ',1)[1]
    if ': ' in v:k,val=v.split(': ',1);cur[k]=val
  if cur:records.append(cur)
 return pd.DataFrame(records)
for acc,name in [('GSE167523','GSE167523_Raw_gene_counts_matrix.txt.gz'),('GSE126848','GSE126848_Gene_counts_raw.txt.gz')]:
 q=O/acc;q.mkdir();record={'cohort':acc};receipt=None
 for i in range(360):
  rr=json.loads((I/'download_receipts.json').read_text()) if (I/'download_receipts.json').exists() else [];found=[z for z in rr if z['name']==name]
  if found and found[0].get('complete_download'):receipt=found[0];break
  if acc=='GSE126848' and (I/'GSE126848_alternative_receipt.json').exists():receipt=json.loads((I/'GSE126848_alternative_receipt.json').read_text());break
  if found and acc!='GSE126848':receipt=found[0];break
  time.sleep(10)
 if not receipt or not receipt.get('complete_download') or receipt.get('error'):
  record.update(status='download_blocked',detail=receipt);rows.append(record);continue
 try:
  raw=pd.read_csv(I/name,sep='\t',index_col=0);assert raw.columns.is_unique and np.isfinite(raw.to_numpy()).all() and raw.to_numpy().min()>=0;assert np.allclose(raw.to_numpy(),np.round(raw.to_numpy()));assert raw.sum(0).gt(0).all();origrows=len(raw);ids=raw.index.astype(str);mapped=pd.Series(ids.str.split('.').str[0],index=ids).map(h) if ids.str.startswith('ENSG').mean()>.5 else pd.Series(ids,index=ids);valid=mapped.notna()&~mapped.duplicated(keep=False);norm=np.log2(raw.div(raw.sum(0),axis=1)*1e6+1).loc[valid.to_numpy()];norm.index=mapped[valid].values;assert norm.index.is_unique;norm.to_csv(q/'log2CPM_full.tsv.gz',sep='\t');meta=soft_meta(acc);meta.to_csv(q/'official_sample_metadata.tsv',sep='\t',index=False);design=pd.DataFrame(index=norm.columns);design.index.name='expression_id';design['clinical_group']='unmatched'
  if acc=='GSE167523':
   clinical=pd.read_excel(I/'GSE167523_Sample_phenotype_correspondence.xlsx');clinical.index=clinical['Sample ID'].map(lambda z:'NAFLD'+str(int(z)).zfill(2));assert clinical.index.is_unique and set(norm.columns)==set(clinical.index);assert set(clinical['GEO accession ID'])==set(meta.GSM);design['GSM']=clinical.loc[design.index,'GEO accession ID'];design['clinical_group']=clinical.loc[design.index,'NAFL/NASH'];record['alignment_rule']='AuthorSampleID integer encoded as provided expression NAFLD01..98; all98columns and GSMs matched one-to-one'
  else:
   bytitle=meta.set_index('title');direct=design.index.isin(bytitle.index)
   if direct.all():design['GSM']=bytitle.loc[design.index,'GSM'].values;design['clinical_group']=bytitle.loc[design.index,'disease'].values;record['alignment_rule']='Exactexpressioncolumn to GEO title'
   elif design.index.isin(meta.GSM).all():mm=meta.set_index('GSM');design['GSM']=design.index;design['clinical_group']=mm.loc[design.index,'disease'].values;record['alignment_rule']='ExactGSM'
   else:record['alignment_rule']='No exactclinicaljoin; overallcolumn-level correlations only; no guessed pairing'
  design.to_csv(q/'sample_design.tsv',sep='\t');r=f".libPaths(c('{P}/R_library','{D}/irGSEA_dependencies_v1/R_project_library_v1','{D}/trajectory_dependencies_v1/R_project_library_v1',.libPaths()));library(GSVA);library(BiocParallel);g<-strsplit(readLines('{D}/F6E_GSE164760_fibrosis_score_v1/HP_HEPATIC_FIBROSIS.gmt'),'\\t');s<-lapply(g,function(z)z[-c(1,2)]);names(s)<-sapply(g,function(z)z[1]);x<-as.matrix(read.delim(gzfile('{q}/log2CPM_full.tsv.gz'),row.names=1,check.names=FALSE));z<-gsva(ssgseaParam(x,s,minSize=1,maxSize=Inf,normalize=TRUE),BPPARAM=SerialParam());write.table(z,'{q}/fibrosis_ssGSEA.tsv',sep='\\t',quote=FALSE,col.names=NA)";(q/'run.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
  with open(q/'GSVA.log','w') as f:proc=subprocess.run(['/usr/bin/Rscript','--vanilla',str(q/'run.R')],env=env,stdout=f,stderr=subprocess.STDOUT)
  assert proc.returncode==0;score=pd.read_csv(q/'fibrosis_ssGSEA.tsv',sep='\t',index_col=0).iloc[0];assert score.index.equals(norm.columns)
  groups=[('all',design.index)]+[(str(k),z.index) for k,z in design.groupby('clinical_group') if k!='unmatched' and len(z)>=5]
  for stratum,ids2 in groups:
   for gene in tfs+[ligand]:
    if gene not in norm.index:continue
    for target in ['HP_HEPATIC_FIBROSIS',ligand,'LDHB','D2HGDH']:
     if target==gene or (target!='HP_HEPATIC_FIBROSIS' and target not in norm.index):continue
     y=score.loc[ids2] if target=='HP_HEPATIC_FIBROSIS' else norm.loc[target,ids2];rho,p=spearmanr(norm.loc[gene,ids2],y);corr.append({'cohort':acc,'stratum':stratum,'gene':gene,'target':target,'samples':len(ids2),'rho':rho,'p':p})
  record.update(status='completed_correlations_only_no_F_labels',raw_genes=origrows,samples=raw.shape[1],unique_mapped_genes=len(norm),excluded_ambiguous_or_unmapped_genes=origrows-len(norm),matched_TFs=len(set(tfs)&set(norm.index)),clinical_group_counts=design.clinical_group.value_counts().to_dict(),F_labels_available=False,complete_input_SHA256=receipt['sha256']);rows.append(record);print(json.dumps(record),flush=True)
 except Exception as e:record.update(status='failed_parse_or_analysis',error=repr(e));rows.append(record);print(json.dumps(record),flush=True)
 pd.DataFrame(rows).to_csv(O/'cohort_coverage.tsv',sep='\t',index=False)
c=pd.DataFrame(corr)
if len(c):
 c['BH_q']=np.nan
 for _,z in c.groupby(['cohort','stratum','target']):
  good=z.p.notna();c.loc[z.index[good],'BH_q']=multipletests(z.loc[good,'p'],method='fdr_bh')[1]
 c.to_csv(O/'extra_TF_ligand_fibrosis_correlations.tsv',sep='\t',index=False)
v={'status':'finished_with_explicit_cohort_status','candidate_sha256':hashlib.sha256(tfpath.read_bytes()).hexdigest(),'cohorts':rows,'ligand':ligand,'scope':'New publicRNAseq counts matched to availableclinicalmetadata. HPOssGSEA and currentTF/PDGFC associations; no F3/F4 AUC or modelrefit without actualFlabels. General NAFL/NASH cohort associations not stage-specific independent clinicalvalidation.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
