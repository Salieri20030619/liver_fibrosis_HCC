from pathlib import Path
import json,collections,numpy as np,pandas as pd,anndata as ad
from scipy.stats import mannwhitneyu,hypergeom
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'continuation_expanded_DE_KEGG_v1';O.mkdir(exist_ok=False);a=ad.read_h5ad(P/'continuation_annotation_evidence_v1/core55_record_counts.h5ad');m=pd.read_csv(P/'continuation_annotation_refined_v1/cell_metadata_reviewed.tsv.gz',sep='\t',index_col=0);assert a.obs_names.equals(m.index)
ks=[q for q in D.glob('*KEGG*') if all((q/f).exists() for f in ['genes.tsv','links.tsv','pathways.tsv'])];assert ks,'No cached complete KEGG official resource';K=sorted(ks)[0];syms=collections.defaultdict(set)
for line in (K/'genes.tsv').read_text().splitlines():
 z=line.split('\t');syms[z[-1].split(';')[0].split(',')[0].strip()].add(z[0])
mapping={k:next(iter(v)) for k,v in syms.items() if len(v)==1};paths=collections.defaultdict(set)
for line in (K/'links.tsv').read_text().splitlines():
 gene,path=line.split('\t');paths[path.replace('path:','')].add(gene)
names=dict(line.split('\t',1) for line in (K/'pathways.tsv').read_text().splitlines());results=[];orows=[];idx={k:i for i,k in enumerate(m.index)}
sub=m[m.reviewed_broad_lineage.isin(['Stromal','Mural','Myeloid'])];scopes={'F3_mesenchymal':sub[sub.reviewed_broad_lineage.isin(['Stromal','Mural'])].index,'F5_myeloid':sub[sub.reviewed_broad_lineage.eq('Myeloid')].index};scopes.update({str(k):g.index for k,g in sub.groupby('refined_cluster')})
for scope,ids in scopes.items():
 q=O/scope;q.mkdir();z=m.loc[ids];raw=[];meta=[]
 for rec,g in z.groupby('donor'):
  if len(g)<10:continue
  assert g.F.nunique()==g.study.nunique()==1;raw.append(np.asarray(a.X[[idx[k] for k in g.index]].sum(0)).ravel());meta.append({'record':rec,'F':int(g.F.iloc[0]),'study':g.study.iloc[0],'cells':len(g)})
 if not raw:continue
 raw=np.vstack(raw);md=pd.DataFrame(meta);md.to_csv(q/'record_metadata.tsv',sep='\t',index=False);cpm=raw/np.maximum(raw.sum(1,keepdims=True),1)*1e6;log=np.log2(cpm+1)
 for study in ['pooled']+sorted(md.study.unique().tolist()):
  sel=np.ones(len(md),bool) if study=='pooled' else md.study.eq(study).to_numpy();high=sel&md.F.ge(3).to_numpy();low=sel&md.F.lt(3).to_numpy();row={'scope':scope,'study':study,'high_records':int(high.sum()),'low_records':int(low.sum())}
  if min(high.sum(),low.sum())<3:row['status']='not_estimable_fewer3_pergroup';results.append(row);continue
  keep=(cpm[sel]>1).sum(0)>=max(3,int(np.ceil(sel.sum()*.1)));test=mannwhitneyu(log[high][:,keep],log[low][:,keep],axis=0,alternative='two-sided',method='asymptotic');d=pd.DataFrame({'gene':a.var_names[keep],'mean_log2CPM_difference':log[high][:,keep].mean(0)-log[low][:,keep].mean(0),'p':test.pvalue});d['BH_q']=multipletests(d.p,method='fdr_bh')[1];d.to_csv(q/(study+'_DE.tsv.gz'),sep='\t',index=False);d['KEGG_ID']=d.gene.map(mapping);dd=d[d.KEGG_ID.notna()].copy();dd=dd[~dd.KEGG_ID.duplicated(keep=False)];universe=set(dd.KEGG_ID);dd[['gene','KEGG_ID']].to_csv(q/(study+'_KEGG_background.tsv'),sep='\t',index=False)
  for direction in ['up','down']:
   sig=dd[(dd.BH_q<.05)&((dd.mean_log2CPM_difference>.25) if direction=='up' else (dd.mean_log2CPM_difference<-.25))];foreground=set(sig.KEGG_ID);rr=[]
   for term,genes in paths.items():
    bg=genes&universe
    if not bg:continue
    hits=bg&foreground;rr.append({'pathway':term,'name':names.get(term,term),'background_genes':len(universe),'pathway_background':len(bg),'foreground_genes':len(foreground),'overlap':len(hits),'p':float(hypergeom.sf(len(hits)-1,len(universe),len(bg),len(foreground))) if foreground else 1.,'overlap_IDs':';'.join(sorted(hits))})
   t=pd.DataFrame(rr);t['BH_q']=multipletests(t.p,method='fdr_bh')[1];t.to_csv(q/(study+'_'+direction+'_KEGG_ORA.tsv'),sep='\t',index=False);orows.append({**row,'direction':direction,'mapped_background':len(universe),'foreground_genes':len(foreground),'significant_pathways':int((t.BH_q<.05).sum())})
  row.update(status='completed',genes=len(d),BH_significant_genes=int((d.BH_q<.05).sum()));results.append(row)
 pd.DataFrame(results).to_csv(O/'DE_comparison_status.tsv',sep='\t',index=False);pd.DataFrame(orows).to_csv(O/'ORA_summary.tsv',sep='\t',index=False);print(scope,'completed',flush=True)
v={'status':'completed','comparisons':len(results),'estimable':sum(r['status']=='completed' for r in results),'KEGG_source':str(K),'design':'Sum rawcounts within final supported record/state or compartment, >=10cells; separate studies plus descriptive pooled; F3/F4 versusF0/1/2 with>=3records/group. CPM>1 in>=max(3,10percent) includedrecords. Two-sided MannWhitney, BH percontrast. ORA foreground BH<.05 and |meanlog2CPMdifference|>.25; tested mappedgenesbackground; up/downseparate; BH perpathwayfamily. Emptyforeground retained asnegative, no nominalp fallback.','limits':'Crossstudy independent people notconfirmed; pooled sourceconfounding; state-level contrasts reflectcapture and composition; no longitudinal/causalclaim. Gene difference is log2CPM mean difference, not count-based foldchange. Noformalplots.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
