import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,datetime,traceback,numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
P=Path(__file__).parent;O=P/'donor_pseudobulk_preparation_v1';O.mkdir(exist_ok=False)
for g in ['GSE244832','GSE136103']:
 try:
  d=O/g;d.mkdir();obs=pd.read_csv(P/'coarse_celltype_review_v2'/(g+'_reviewed_cell_metadata.tsv.gz'),sep='\t',index_col=0,keep_default_na=False);a=ad.read_h5ad(P/g/'qc_pass_counts.h5ad',backed='r');assert obs.index.tolist()==a.obs_names.tolist()
  fields=['donor','library','sorting','coarse_label_reviewed'];keys=obs[fields].astype(str).agg('|'.join,axis=1);codes,levels=pd.factorize(keys,sort=True);sums=np.zeros((len(levels),a.n_vars),dtype='float64');cells=np.bincount(codes,minlength=len(levels))
  for start in range(0,a.n_obs,5000):
   end=min(start+5000,a.n_obs);xx=a.X[start:end,:];inc=sp.coo_matrix((np.ones(end-start),(codes[start:end],np.arange(end-start))),shape=(len(levels),end-start)).tocsr();sums+=(inc@xx).toarray()
  assert np.isfinite(sums).all() and (sums>=0).all() and np.allclose(sums,np.rint(sums));meta=[]
  for i,key in enumerate(levels):
   group=obs.iloc[np.flatnonzero(codes==i)];row={'pseudobulk_id':g+'_'+str(i+1),'cells':int(cells[i]),'library_sum':float(sums[i].sum()),'statistical_unit':'donor; libraries and sorting strata are not independent donors'}
   for c in fields+['coarse_review_status','condition_conflict','disease_comparison_eligible','condition_metadata_original','condition_Excel_original','F_stage_from_Figure4C','disease status','cause of liver disease']:
    if c in group:row[c]=';'.join(sorted(group[c].astype(str).unique()))
   meta.append(row)
  m=pd.DataFrame(meta).set_index('pseudobulk_id');var=a.var.copy();a.file.close();out=ad.AnnData(sp.csr_matrix(sums),obs=m,var=var);out.uns['scope']='Raw library-by-coarse-label sums; no differential test; preserve sorting; aggregate repeated libraries by donor only within compatible design';out.write_h5ad(d/'library_stratified_pseudobulk_counts.h5ad',compression='gzip');m.to_csv(d/'pseudobulk_sample_manifest.tsv',sep='\t')
  result={'input_cells':len(obs),'pseudobulk_rows':len(m),'donors':obs.donor.nunique(),'sum_cell_counts_matches':int(cells.sum())==len(obs),'integer_nonnegative_counts':True,'genes':len(var),'differential_model_run':False,'unresolved_labels_retained_as_separate_groups':True,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};(d/'result.json').write_text(json.dumps(result));print(json.dumps({'cohort':g,**result}),flush=True)
 except Exception:(O/(g+'_error.txt')).write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
