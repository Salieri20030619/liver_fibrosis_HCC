import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad,datetime,traceback
P=Path(__file__).parent;O=P/'expanded_marker_review_v1';O.mkdir(exist_ok=False)
markers='PTPRC TRAC CD3D CD3E CD247 LCK CD8A CD8B NKG7 GNLY KLRD1 FCGR3A NCAM1 KLRB1 TRDC MS4A1 CD79A CD79B MZB1 JCHAIN XBP1 SDC1 LYZ CD68 CD163 MARCO C1QA C1QB CSF1R FCER1A CLEC9A CD1C IL3RA GZMB KRT8 KRT18 KRT7 KRT19 EPCAM SOX9 HNF4A TTR APOA1 ASGR1 ALB PDGFRB PDGFRA DES RBP1 LRAT COL1A1 COL3A1 DCN RGS5 MYH11 ACTA2 PECAM1 VWF KDR EMCN CLEC4G STAB2 PLVAP ESAM'.split()
for g in ['GSE244832','GSE136103']:
 try:
  a=ad.read_h5ad(P/g/'qc_pass_counts.h5ad',backed='r');b=ad.read_h5ad(P/'F2_atlas_preparation_v2'/g/'embedding_checkpoint.h5ad',backed='r');assert a.obs_names.equals(b.obs_names);obs=b.obs.copy();b.file.close();syms=a.var.gene_symbol_original.astype(str).to_numpy();coverage=[{'marker':m,'matching_features':int((syms==m).sum())} for m in markers];pd.DataFrame(coverage).to_csv(O/(g+'_marker_coverage.tsv'),sep='\t',index=False);valid=[m for m in markers if (syms==m).sum()==1];idx=[int(np.flatnonzero(syms==m)[0]) for m in valid];val=np.empty((a.n_obs,len(idx)),dtype='float32');det=np.empty(val.shape,dtype=bool)
  for start in range(0,a.n_obs,5000):
   end=min(start+5000,a.n_obs);xx=a.X[start:end,:];total=np.asarray(xx.sum(axis=1)).ravel();sel=xx[:,idx].toarray();val[start:end]=np.log1p(sel/np.maximum(total[:,None],1)*1e4);det[start:end]=sel>0
  a.file.close();rows=[];joint=[]
  for group,ix in obs.groupby(['cluster_preliminary','library'],observed=True).indices.items():
   for j,m in enumerate(valid):rows.append({'cluster':group[0],'library':group[1],'cells':len(ix),'marker':m,'mean_log1p_cp10k':float(val[ix,j].mean()),'fraction_detected':float(det[ix,j].mean())})
   for m,n in [('ALB','ASGR1'),('ALB','EPCAM'),('CD3D','NKG7'),('CD3D','MS4A1'),('PTPRC','EPCAM')]:
    if m in valid and n in valid:joint.append({'cluster':group[0],'library':group[1],'cells':len(ix),'pair':m+' & '+n,'co_detected_fraction':float((det[ix,valid.index(m)]&det[ix,valid.index(n)]).mean())})
  pd.DataFrame(rows).to_csv(O/(g+'_expanded_markers_by_cluster_library.tsv.gz'),sep='\t',index=False);pd.DataFrame(joint).to_csv(O/(g+'_co_detection_review.tsv'),sep='\t',index=False);(O/(g+'_result.json')).write_text(json.dumps({'cells':len(obs),'markers_unique':len(valid),'rows':len(rows),'no_automatic_identity_assignment':True,'co_detection_not_doublet_diagnosis':True,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}));print(g+' complete',flush=True)
 except Exception:(O/(g+'_error.txt')).write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
