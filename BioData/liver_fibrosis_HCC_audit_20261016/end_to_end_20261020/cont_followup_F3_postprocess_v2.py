from pathlib import Path
import time,json,numpy as np,pandas as pd,collections
from scipy.io import mmread
from scipy.stats import spearmanr,hypergeom
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';T=P/'continuation_F3_CytoTRACE_serial_resume_v2';I=P/'continuation_followup_inputs_v1';O=P/'continuation_F3_trajectory_checks_v2';O.mkdir(exist_ok=False)
for i in range(1080):
 if (T/'validation.json').exists():break
 time.sleep(10)
else:raise RuntimeError('Trajectory job has not reached a terminal validation within3h')
tv=json.loads((T/'validation.json').read_text());rows=[];associations=[]
K=D/'F3F4_KEGG_function_v2';syms=collections.defaultdict(set)
for line in (K/'genes.tsv').read_text().splitlines():
 z=line.split('\t');syms[z[-1].split(';')[0].split(',')[0].strip()].add(z[0])
mapping={k:next(iter(v)) for k,v in syms.items() if len(v)==1};paths=collections.defaultdict(set)
for line in (K/'links.tsv').read_text().splitlines():
 g,p=line.split('\t');paths[p.replace('path:','')].add(g)
names=dict(line.split('\t',1) for line in (K/'pathways.tsv').read_text().splitlines())
for branch in ['F3']:
 q=T;out=O/branch;out.mkdir()
 if not (q/'COMPLETE.txt').exists():rows.append({'branch':branch,'status':'upstream_failed','error':tv.get('error_tail')});continue
 m=pd.read_csv(q/'ordered_cell_metadata.tsv',sep='\t',index_col=0);base=pd.read_csv(I/branch/'metadata.tsv',sep='\t',index_col=0);assert m.index.equals(base.index);assert (m.refined_cluster==base.refined_cluster).all();assert np.isfinite(m[['Pseudotime','CytoTRACE']]).all().all();de=pd.read_csv(q/'pseudotime_dynamic_genes.tsv',sep='\t',index_col=0);assert len(de)==3000;de.to_csv(out/'dynamic_gene_results.tsv.gz',sep='\t');means=m.groupby(['study','donor','F']).agg(cells=('Pseudotime','size'),pseudotime_mean=('Pseudotime','mean'),CytoTRACE_mean=('CytoTRACE','mean')).reset_index();means.to_csv(out/'record_means.tsv',sep='\t',index=False)
 for study in ['pooled']+sorted(means.study.unique().tolist()):
  z=means if study=='pooled' else means[means.study.eq(study)]
  for val in ['pseudotime_mean','CytoTRACE_mean']:
   r,p=spearmanr(z.F,z[val]);associations.append({'branch':branch,'study':study,'metric':val,'records':len(z),'rho_F':r,'p':p})
 pd.crosstab(m.State,m.refined_cluster).to_csv(out/'trajectory_state_annotation_counts.tsv',sep='\t');pd.crosstab(m.State,m.study).to_csv(out/'trajectory_state_study_counts.tsv',sep='\t');pd.crosstab(m.State,m.F).to_csv(out/'trajectory_state_F_counts.tsv',sep='\t')
 # Numeric trend matrix only: global pseudotime quantile bins, equal weight to each observed record/bin with >=5cells.
 genes=(I/branch/'genes.txt').read_text().splitlines();x=mmread(str(I/branch/'raw_counts.mtx')).tocsr();total=np.asarray(x.sum(0)).ravel();take=[genes.index(g) for g in de.index];y=x[take].toarray();y=np.log1p(y/np.maximum(total,1)*1e4);bins=pd.qcut(m.Pseudotime,20,labels=False,duplicates='drop');profiles=[];bn=[]
 for b in sorted(bins.unique()):
  sub=m[bins.eq(b)];vals=[]
  for record,g in sub.groupby('donor'):
   if len(g)>=5:vals.append(y[:,m.index.get_indexer(g.index)].mean(1))
  if vals:profiles.append(np.vstack(vals).mean(0));bn.append({'bin':int(b),'records_with5cells':len(vals),'cells_in_bin':len(sub),'pseudotime_median':float(sub.Pseudotime.median())})
 pd.DataFrame(np.vstack(profiles).T,index=de.index,columns=[r['bin'] for r in bn]).to_csv(out/'record_equal_pseudotime_binned_log1pCP10K.tsv.gz',sep='\t');pd.DataFrame(bn).to_csv(out/'bin_coverage.tsv',sep='\t',index=False)
 dd=de.copy();dd['KEGG_ID']=dd.index.map(mapping);dd=dd[dd.KEGG_ID.notna()];dd=dd[~dd.KEGG_ID.duplicated(keep=False)];universe=set(dd.KEGG_ID);fg=set(dd.loc[dd.qval.lt(.05),'KEGG_ID']);rr=[]
 for term,gs in paths.items():
  bg=gs&universe
  if not bg:continue
  hit=bg&fg;rr.append({'pathway':term,'name':names.get(term,term),'universe':len(universe),'foreground':len(fg),'pathway_background':len(bg),'overlap':len(hit),'p':float(hypergeom.sf(len(hit)-1,len(universe),len(bg),len(fg))) if fg else 1,'overlap_IDs':';'.join(sorted(hit))})
 rr=pd.DataFrame(rr);rr['BH_q']=multipletests(rr.p,method='fdr_bh')[1];rr.to_csv(out/'dynamic_gene_KEGG_ORA.tsv',sep='\t',index=False);rows.append({'branch':branch,'status':'completed','cells':len(m),'source_records':m.donor.nunique(),'reviewed_states':m.refined_cluster.nunique(),'trajectory_states':m.State.nunique(),'tested_genes':len(de),'dynamic_q05':int(de.qval.lt(.05).sum()),'dynamic_test_status_counts':de.status.value_counts().to_dict() if 'status' in de else {},'dynamic_KEGG_q05':int(rr.BH_q.lt(.05).sum()),'all_cell_IDs_labels_match_final_inputs':True})
a=pd.DataFrame(associations)
if len(a):a['BH_q']=multipletests(a.p,method='fdr_bh')[1];a.to_csv(O/'record_F_associations.tsv',sep='\t',index=False)
v={'status':'completed' if all(r['status']=='completed' for r in rows) else 'partial_upstream_failure','branches':rows,'limits':'Dynamic tests and ORA reflect cell-level pseudotime associations, not independent-record progression. Record-equal binned profiles are a numeric descriptive reconstruction, not original unspecified smoothing. Per-study F associations reported; no trajectory direction establishes causality or HCC lineage.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
