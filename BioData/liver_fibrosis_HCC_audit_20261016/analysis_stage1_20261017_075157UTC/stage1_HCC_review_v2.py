import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,time,itertools,traceback,numpy as np,pandas as pd,anndata as ad,scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;A=P/'F2_HCC_GSE189175_v2';O=P/'HCC_marker_pair_review_v2';O.mkdir(exist_ok=False)
try:
 for _ in range(240):
  state=json.loads((A/'state.json').read_text())
  if state['status']=='complete':break
  if state['status']=='failed':raise RuntimeError(state['error'])
  time.sleep(15)
 else:raise TimeoutError('Atlas not complete after 1 hour')
 b=ad.read_h5ad(A/'GSE189175/atlas_HVG_preparation.h5ad',backed='r');obs=b.obs.copy();b.file.close();path=json.loads((P/'downstream_QC_input_overrides_v2.json').read_text())['GSE189175'];a=ad.read_h5ad(path);assert a.obs_names.equals(obs.index);assert obs.donor.nunique()==3 and len(obs)==31410;assert not obs.mt_metric_available.astype(bool).any();X=a.X.tocsr();symbols=a.var.gene_symbol_original.astype(str);markers='ALB ASGR1 CYP3A4 APOA1 TTR KRT8 KRT18 GPC3 AFP EPCAM KRT7 KRT19 SOX9 PTPRC LYZ LST1 TYROBP AIF1 CD68 CD163 MARCO C1QA C1QB S100A8 S100A9 FCN1 VCAN CD14 FCGR3A CLEC9A CD1C FCER1A CLEC4C GZMB CD3D CD3E TRAC IL7R CD8A NKG7 GNLY KLRD1 MS4A1 CD79A CD79B MZB1 JCHAIN PECAM1 VWF KDR EMCN STAB2 LYVE1 PROX1 COL1A1 COL1A2 COL3A1 DCN LUM PDGFRA RGS5 CSPG4 ACTA2 TAGLN MYH11 DES MKI67 TOP2A'.split();valid=[];missing=[];values=[];total=np.asarray(X.sum(axis=1)).ravel()
 for m in markers:
  ix=np.flatnonzero(symbols.to_numpy()==m)
  if not len(ix):missing.append(m);continue
  v=np.asarray(X[:,ix].sum(axis=1)).ravel();values.append(v);valid.append(m)
 raw=np.array(values).T;norm=np.log1p(raw/total[:,None]*1e4);clusters=sorted(obs.cluster_preliminary.astype(str).unique(),key=int);rows=[]
 for cl in clusters:
  ix=obs.cluster_preliminary.astype(str).eq(cl).to_numpy()
  for j,m in enumerate(valid):rows.append({'cluster':cl,'marker':m,'cells':int(ix.sum()),'mean_log1p_CP10k':float(norm[ix,j].mean()),'fraction_detected':float((raw[ix,j]>0).mean())})
 t=pd.DataFrame(rows);t.to_csv(O/'cluster_marker_evidence.tsv',sep='\t',index=False);mean=t.pivot(index='cluster',columns='marker',values='mean_log1p_CP10k').reindex(index=clusters,columns=valid);frac=t.pivot(index='cluster',columns='marker',values='fraction_detected').reindex(index=clusters,columns=valid);z=(mean-mean.mean())/mean.std(ddof=0).replace(0,1);fig,ax=plt.subplots(figsize=(22,max(5,len(clusters)*.35)));xx,yy=np.meshgrid(np.arange(len(valid)),np.arange(len(clusters)));im=ax.scatter(xx.ravel(),yy.ravel(),s=frac.to_numpy().ravel()*75,c=z.to_numpy().ravel(),cmap='RdBu_r',vmin=-2,vmax=2);ax.set_xticks(range(len(valid)),valid,rotation=90,fontsize=8);ax.set_yticks(range(len(clusters)),clusters);ax.set_ylabel('Unassigned cluster');ax.set_title('GSE189175 marker evidence | 3 paired donors | dot area = detection; color = marker z-score');fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(O/'HCC_cluster_marker_dotplot.png',dpi=160);plt.close(fig)
 sample=obs[['donor','library','tissue','Condition']].drop_duplicates();assert sample.library.nunique()==6 and sample.groupby('donor').tissue.nunique().eq(2).all();sample.to_csv(O/'sample_pairing.tsv',sep='\t',index=False);counts=pd.crosstab(obs.library,obs.cluster_preliminary.astype(str)).reindex(columns=clusters,fill_value=0);counts.to_csv(O/'cluster_by_library_counts.tsv',sep='\t');fractions=counts.div(counts.sum(axis=1),axis=0);fractions.to_csv(O/'cluster_by_library_fractions.tsv',sep='\t');tests=[];fig,axs=plt.subplots(int(np.ceil(len(clusters)/4)),4,figsize=(15,int(np.ceil(len(clusters)/4))*3),squeeze=False);signs=np.array(list(itertools.product([-1,1],repeat=3)))
 for j,cl in enumerate(clusters):
  dif=[];ax=axs.flat[j]
  for donor,m in sample.groupby('donor'):
   adjacent=m.loc[m.tissue.eq('AdjNonTumor'),'library'].iloc[0];tumor=m.loc[m.tissue.eq('Tumor'),'library'].iloc[0];v=[fractions.loc[adjacent,cl],fractions.loc[tumor,cl]];dif.append(v[1]-v[0]);ax.plot([0,1],v,'o-',label=str(donor),alpha=.8)
  dif=np.array(dif);effect=dif.mean();p=float((np.abs((signs*dif).mean(axis=1))>=abs(effect)-1e-14).mean());tests.append({'cluster':cl,'independent_pairs':3,'mean_paired_fraction_difference_tumor_minus_adjacent':effect,'exact_two_sided_sign_flip_p':p,'positive_pairs':int((dif>0).sum()),'negative_pairs':int((dif<0).sum())});ax.set_title('Cluster '+cl);ax.set_xticks([0,1],['Adjacent','Tumor']);ax.set_ylabel('Cell fraction')
 for ax in list(axs.flat)[len(clusters):]:ax.axis('off')
 axs.flat[0].legend(title='Patient',fontsize=7);fig.suptitle('GSE189175 paired composition | adjacent tissue is not healthy liver');fig.tight_layout();fig.savefig(O/'HCC_paired_cluster_fractions.png',dpi=160);plt.close(fig);tt=pd.DataFrame(tests);pv=tt.exact_two_sided_sign_flip_p.to_numpy();order=np.argsort(pv);q=np.empty(len(pv));q[order]=np.minimum(1,np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1]);tt['BH_q']=q;tt.to_csv(O/'paired_cluster_fraction_tests.tsv',sep='\t',index=False)
 groups=obs.groupby(['donor','library','tissue','cluster_preliminary'],observed=True).indices;pb=[];md=[]
 for key,ix in groups.items():pb.append(sp.csr_matrix(X[ix].sum(axis=0)));md.append(dict(zip(['donor','library','tissue','cluster'],map(str,key)),cells=len(ix)))
 pb=sp.vstack(pb).tocsr();assert int(pb.sum())==int(X.sum());sp.save_npz(O/'cluster_sample_raw_pseudobulk.npz',pb);pd.DataFrame(md).to_csv(O/'pseudobulk_rows.tsv',sep='\t',index=False);a.var.to_csv(O/'pseudobulk_genes.tsv',sep='\t');r={'cells':len(obs),'genes':a.n_vars,'clusters':len(clusters),'donors':3,'tissues':6,'marker_genes_available':len(valid),'missing_markers':missing,'pairs_verified':True,'pseudobulk_count_conservation':True,'fraction_tests_q_below05':int((q<.05).sum()),'minimum_possible_two_sided_sign_flip_p_with_3_pairs':.25,'MT_QC':'unavailable; corrected NA metric retained; no mitochondrial filtering applied','labels':'Clusters unassigned pending marker review; no malignant identity inferred from tumor origin','source_condition_labels':sample.Condition.unique().tolist()};(O/'validation.json').write_text(json.dumps(r,indent=2));print(json.dumps(r))
except Exception:(O/'error.txt').write_text(traceback.format_exc());raise
