import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,datetime,sys,gzip,tarfile,io,traceback,hashlib
import numpy as np,pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from threadpoolctl import threadpool_limits
R=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');O=Path(__file__).parent
MARKERS={'Hepatocyte':['ALB','ASGR1','CYP3A4'],'Endothelial':['PECAM1','VWF','STAB2'],'HSC_fibroblast':['DCN','COL3A1','COL1A1'],'T_cell':['CD3D','CD3E','CD2'],'B_cell':['MS4A1','CD79A','CD79B'],'Myeloid':['LYZ','CD163','MARCO'],'Cholangiocyte':['KRT7','KRT19','EPCAM'],'NK':['NKG7','GNLY','KLRD1'],'SMC_pericyte':['RGS5','MYH11','ACTA2']}
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def state(task,status,**kw):
 x={'task':task,'status':status,'utc':stamp(),**kw};(O/(task+'.state.json')).write_text(json.dumps(x,indent=2,default=str));print(json.dumps(x,default=str),flush=True)
def qc(g,X,genes,obs):
 import anndata as ad,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 d=O/g;d.mkdir(exist_ok=True);state('singlecell','qc',cohort=g,cells=X.shape[0],genes=X.shape[1])
 X=X.tocsr().astype(np.float32);X.eliminate_zeros();genes=np.asarray(genes,dtype=str)
 obs=obs.copy();obs.index=obs.index.astype(str);obs['study']=g
 obs['total_counts_recomputed']=np.asarray(X.sum(axis=1)).ravel();obs['n_genes_recomputed']=np.diff(X.indptr)
 mt=np.char.startswith(np.char.upper(genes),'MT-');obs['pct_mt_recomputed']=np.asarray(X[:,mt].sum(axis=1)).ravel()/np.maximum(obs.total_counts_recomputed,1)*100
 thresholds=[];obs['qc_fail_3mad']=False
 for lib,ix in obs.groupby('library',sort=False).indices.items():
  for metric,mode in [('total_counts_recomputed','two_sided'),('n_genes_recomputed','two_sided'),('pct_mt_recomputed','upper')]:
   vals=obs.iloc[ix][metric].to_numpy();med=float(np.median(vals));mad=float(np.median(np.abs(vals-med)));low=max(0,med-3*mad);high=med+3*mad
   # Zero MAD cannot define a reliable outlier fence; preserve values and flag for review.
   fail=(vals>high)|((vals<low)&(mode=='two_sided')) if mad>0 else np.zeros(len(vals),bool)
   col='fail_'+metric
   if col not in obs:obs[col]=False
   obs.iloc[ix,obs.columns.get_loc(col)]=fail
   thresholds.append({'library':lib,'metric':metric,'mode':mode,'median':med,'MAD_unscaled':mad,'lower':low,'upper':high,'zero_MAD_no_filter':mad==0,'flagged':int(fail.sum())})
 failcols=[c for c in obs if c.startswith('fail_')];obs['qc_fail_3mad']=obs[failcols].any(axis=1);obs['qc_pass']=~obs.qc_fail_3mad & (obs.total_counts_recomputed>0)
 obs.to_csv(d/'cell_QC_and_source_metadata.tsv.gz',sep='\t');pd.DataFrame(thresholds).to_csv(d/'QC_thresholds_per_library.tsv',sep='\t',index=False)
 summary=obs.groupby(['donor','library'],dropna=False).agg(input_cells=('qc_pass','size'),pass_cells=('qc_pass','sum')).reset_index();summary.to_csv(d/'QC_by_donor_library.tsv',sep='\t',index=False)
 fig,axs=plt.subplots(1,3,figsize=(12,3.5))
 for ax,col in zip(axs,['total_counts_recomputed','n_genes_recomputed','pct_mt_recomputed']):
  v=obs[col];ax.hist(v,bins=80,color='#4b7892');ax.set_title(col);ax.set_yscale('log')
 fig.suptitle(g+' | input QC distributions; 3 MAD thresholds per library');fig.tight_layout();fig.savefig(d/'F2_preparation_QC.png',dpi=160);plt.close(fig)
 keep=obs.qc_pass.to_numpy();Xp=X[keep].copy();obsp=obs.loc[keep].copy();nonzero=np.asarray(Xp.sum(axis=0)).ravel()>0
 var=pd.DataFrame({'gene_symbol_original':genes,'source_feature_index_1based':np.arange(len(genes))+1},index=pd.Index(genes));a=ad.AnnData(Xp[:,nonzero],obs=obsp,var=var.loc[nonzero].copy());a.var_names_make_unique()
 # Marker evidence only: no forced cell-type assignment and no tumor identity inferred from one marker.
 totals=np.asarray(a.X.sum(axis=1)).ravel();norm=sp.diags(1e4/np.maximum(totals,1)).dot(a.X).tocsr();norm.data=np.log1p(norm.data)
 groupcols=['library']+(['seurat_cluster_original'] if 'seurat_cluster_original' in obsp else [])
 evidence=[];available=[]
 syms=a.var.gene_symbol_original.to_numpy()
 for lineage,markers in MARKERS.items():
  for marker in markers:
   jj=np.flatnonzero(syms==marker);available.append({'program':lineage,'marker':marker,'matched_features':len(jj)})
   if len(jj)!=1:continue
   v=norm[:,jj[0]].toarray().ravel();raw=a.X[:,jj[0]].toarray().ravel()
   for group,ix in obsp.groupby(groupcols,observed=True).indices.items():evidence.append({'group':str(group),'program':lineage,'marker':marker,'cells':len(ix),'mean_log1p_cp10k':float(v[ix].mean()),'fraction_detected':float((raw[ix]>0).mean())})
 pd.DataFrame(evidence).to_csv(d/'F2G_marker_evidence_by_source_group.tsv',sep='\t',index=False);pd.DataFrame(available).to_csv(d/'marker_coverage.tsv',sep='\t',index=False)
 del norm,Xp
 # Keep untouched counts in X; only QC removal is applied. Embedding is a separate stage.
 for c in a.obs:
  if a.obs[c].dtype==object:a.obs[c]=a.obs[c].fillna('').astype(str)
 a.uns['analysis_status']='QC complete; marker evidence only; no final cell identity or integrated atlas'
 a.uns['qc_reconstruction']='Baseline 3MAD; library-wise raw metrics; two-sided counts/genes, upper mitochondrial; grouping/scale/tail not specified in baseline, reconstructed explicitly; zero MAD flagged without filter'
 a.write_h5ad(d/'qc_pass_counts.h5ad',compression='gzip')
 result={'cohort':g,'input_cells':len(obs),'qc_pass_cells':int(keep.sum()),'donors':int(obs.donor.nunique()),'genes_input':len(genes),'expressed_genes_retained':int(nonzero.sum()),'mt_genes':int(mt.sum()),'status':'QC_and_marker_evidence_complete','identity_status':'not_final','outputs':['qc_pass_counts.h5ad','cell_QC_and_source_metadata.tsv.gz','QC_thresholds_per_library.tsv','F2G_marker_evidence_by_source_group.tsv']}
 (d/'result.json').write_text(json.dumps(result,indent=2));state('singlecell','cohort_complete',**{k:v for k,v in result.items() if k!='status'})
def singlecell():
 tasks=json.loads((O/'analysis_manifest.json').read_text())['singlecell']
 for t in tasks:
  g=t['cohort'];d=O/g
  if (d/'result.json').exists():continue
  try:
   state('singlecell','loading_counts',cohort=g,path=t['matrix'])
   if g=='GSE244832':
    base=R/'round2_processed_20261017';inv=json.loads((base/'complete_small_member_field_audit.json').read_text());paths={Path(x['member']).name:Path(x['path']) for x in inv}
    genes=pd.read_csv(paths['hLIVER_genes.csv'],header=None)[0].astype(str).to_numpy();obs=pd.read_csv(base/'cell_donor_audit/GSE244832_cell_matrix_order_donor_F.tsv',sep='\t',keep_default_na=False).set_index('cell_ID');obs['donor']=obs.Donor_ID;obs['library']=obs['orig.ident_original'];obs['tissue']='liver';obs['sorting']='not_explicit';obs['disease_comparison_eligible']=~obs.condition_conflict.astype(str).str.lower().eq('true')
    with tarfile.open(t['matrix'],'r:gz') as tf:
     m=next(m for m in tf.getmembers() if m.name.endswith('hLIVER_counts.mtx.gz'));X=mmread(io.BytesIO(tf.extractfile(m).read())).T.tocsr()
    qc(g,X,genes,obs)
   elif g in ['GSE136103','GSE212837']:
    meta=pd.read_csv(R/(g+'_library_donor_mapping.tsv'),sep='\t',keep_default_na=False).set_index('GSM');xs=[];ob=[];genes=None
    with tarfile.open(t['matrix'],'r:') as tf:
     members={m.name:m for m in tf.getmembers() if m.isfile()}
     for name,m in members.items():
      if not name.endswith('_matrix.mtx.gz'):continue
      gsm=Path(name).name.split('_')[0];s=meta.loc[gsm]
      if str(s.include_human_liver).lower()!='true':continue
      pre=name[:-len('_matrix.mtx.gz')];gn=next(n for n in [pre+'_genes.tsv.gz',pre+'_features.tsv.gz'] if n in members);bn=pre+'_barcodes.tsv.gz'
      feats=[x.split('\t') for x in gzip.decompress(tf.extractfile(gn).read()).decode().splitlines()];gns=np.array([x[1] if len(x)>1 else x[0] for x in feats]);bar=gzip.decompress(tf.extractfile(bn).read()).decode().splitlines()
      if genes is None:genes=gns
      else:assert np.array_equal(genes,gns),'feature order differs; stop instead of positional merge'
      xx=mmread(io.BytesIO(gzip.decompress(tf.extractfile(m).read()))).T.tocsr();oo=pd.DataFrame(index=[gsm+'|'+b for b in bar])
      for k,v in s.items():oo[k]=str(v)
      oo['GSM']=gsm;oo['donor']=str(s.derived_donor_id);oo['library']=gsm;oo['tissue']='liver';oo['sorting']=s.get('cell subtype','not_explicit');oo['disease_comparison_eligible']=True
      xs.append(xx);ob.append(oo);state('singlecell','library_loaded',cohort=g,GSM=gsm,cells=len(bar))
    X=sp.vstack(xs,format='csr');del xs;qc(g,X,genes,pd.concat(ob))
   elif g=='GSE189175':
    genes=[x.split('\t')[1] for x in gzip.open(R/'GSE189175_counts_features.tsv.gz','rt').read().splitlines()];bar=gzip.open(R/'GSE189175_counts_barcodes.tsv.gz','rt').read().splitlines();obs=pd.read_csv(R/'GSE189175_meta_data.tsv.gz',sep='\t',keep_default_na=False,index_col=0).loc[bar];obs['donor']=obs.Patient.astype(str);obs['library']=obs.SampleID.astype(str);obs['sorting']='not_explicit';obs['disease_comparison_eligible']=True
    X=mmread(gzip.open(t['matrix'],'rb')).T.tocsr();qc(g,X,genes,obs)
   del X
  except Exception:
   d.mkdir(exist_ok=True);(d/'error.txt').write_text(traceback.format_exc());state('singlecell','cohort_failed',cohort=g,error=traceback.format_exc())
 state('singlecell','finished_queue')
def bulk():
 from sklearn.decomposition import PCA
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 tasks=json.loads((O/'analysis_manifest.json').read_text())['bulk']
 for t in tasks:
  g=t['cohort'];d=O/g;d.mkdir(exist_ok=True)
  try:
   state('bulk','loading',cohort=g)
   a=pd.read_csv(t['matrix'],sep=t['sep'],index_col=0);m=pd.read_csv(t['metadata'],sep='\t',keep_default_na=False);key=t['join_key'];m=m.set_index(key).loc[a.columns];m.to_csv(d/'sample_manifest_expression_order.tsv',sep='\t')
   a.index=a.index.astype(str);gene=pd.DataFrame({'source_gene_id':a.index,'ensembl_id_without_version':a.index.str.replace(r'\.\d+$','',regex=True),'symbol_status':'not_invented; annotation_pending_or_source_symbol'});gene.to_csv(d/'gene_identifier_inventory.tsv',sep='\t',index=False)
   cp=np.log1p(a.to_numpy(dtype=float)/np.maximum(a.sum(axis=0).to_numpy(),1)[None,:]*1e6);v=np.var(cp,axis=1);sel=np.argsort(v)[-min(3000,len(v)):];pca=PCA(n_components=2,random_state=0).fit_transform(cp[sel].T);pd.DataFrame(pca,index=a.columns,columns=['PC1','PC2']).to_csv(d/'bulk_QC_PCA.tsv',sep='\t')
   fig,ax=plt.subplots(figsize=(6,5));ax.scatter(pca[:,0],pca[:,1],s=18);ax.set(xlabel='PC1',ylabel='PC2',title=g+' | log1p CPM PCA: preparation, not classifier');fig.tight_layout();fig.savefig(d/'bulk_preparation_PCA.png',dpi=160);plt.close(fig)
   pd.DataFrame({'sample':a.columns,'library_sum':a.sum(axis=0).to_numpy(),'detected_genes':(a>0).sum(axis=0).to_numpy()}).to_csv(d/'bulk_sample_QC.tsv',sep='\t',index=False)
   x={'genes':len(a),'samples':a.shape[1],'metadata_order_joined':True,'normalization_for_QC_only':'log1p CPM','fractional_expected_counts_preserved':t.get('fractional_expected_counts',False),'clinical_model_run':False};(d/'result.json').write_text(json.dumps(x,indent=2));state('bulk','cohort_complete',cohort=g,**x)
  except Exception:state('bulk','cohort_failed',cohort=g,error=traceback.format_exc());(d/'error.txt').write_text(traceback.format_exc())
 state('bulk','finished_queue')
if __name__=='__main__':
 with threadpool_limits(limits=4):
  {'singlecell':singlecell,'bulk':bulk}[sys.argv[1]]()
