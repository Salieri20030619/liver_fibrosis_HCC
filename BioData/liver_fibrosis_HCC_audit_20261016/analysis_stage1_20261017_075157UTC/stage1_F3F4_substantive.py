import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
from pathlib import Path
import json,datetime,traceback,subprocess,shutil,importlib.util
import pandas as pd,numpy as np,anndata as ad,scipy.sparse as sp
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'F3F4_donor_analysis_v1';O.mkdir(exist_ok=False);S=P/'F2_Roe_F3_F4_subclusters_v1'
def record(task,**kw):
 v={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'task':task,**kw};(O/'state.json').write_text(json.dumps(v));print(json.dumps(v),flush=True)
def bh(p):
 n=len(p);order=np.argsort(p);q=np.empty(n);q[order]=np.minimum(1,np.minimum.accumulate((np.asarray(p)[order]*n/np.arange(1,n+1))[::-1])[::-1]);return q
manual={
'GSE244832_Mesenchymal':{'0':'ECM_stromal','1':'ECM_stromal','2':'ECM_stromal','3':'Stromal_mural_mixed','4':'Stromal_unresolved','5':'ECM_stromal','6':'Stromal_unresolved','7':'Unresolved_stromal_myeloid_markers','8':'Contractile_mural'},
'GSE136103_Mesenchymal':{str(i):('ECM_stromal' if i in [4,9] else 'Unresolved_mixed' if i in [10,14,15,17] else 'Stromal_mural_mixed' if i in [1,8] else 'Contractile_mural') for i in range(18)},
'GSE244832_Myeloid':{str(i):('Unresolved_low_marker' if i==8 else 'Macrophage_like') for i in range(9)},
'GSE136103_Myeloid':{str(i):('Macrophage_like' if i in [1,3,9] else 'Unresolved_mixed' if i in [8,13] else 'Myeloid_unresolved' if i==11 else 'Inflammatory_myeloid') for i in range(14)}}
review=[]
for task,groups in manual.items():
 for cl,label in groups.items():review.append({'task':task,'cluster':cl,'broad_review_label':label,'basis':'Manual joint review of ECM/contractile/HSC and myeloid/lymphoid/epithelial marker detection plus donor distribution; no numerical assignment threshold','not_inferred':'No DTNA subtype, no RUNX2 axis, no quiescent-HSC root, no doublet diagnosis'})
pd.DataFrame(review).to_csv(O/'manual_broad_identity_review.tsv',sep='\t',index=False)
manifest={'panel_mapping':{'F3A/B_F4A/B':'existing UMAP/marker plots with finite manual broad review','F3C/D_F4C':'within-study within-sorting subcluster Ro/e and donor fractions','F3E':'donor-aggregated compartment differential expression; broad mesenchymal not asserted pure HSC','F3F_F4D-G':'KEGG enrichment of donor DE in next step','F3G-J_F5':'CytoTRACE/Monocle2 pending dependency and biologically supported root'},'contrasts':'GSE244832 F3+F4 versus F0, D9/D10 conflicts excluded; GSE136103 fibrotic versus healthy separately by original sorting','statistics':'Sum raw counts over each donor x sorting x compartment, combine library repeats; log2(CPM+1), two-sided Mann-Whitney rank-sum with tie correction (asymptotic); all source genes tested, no new expression or cell-count filter; BH reported alongside raw p; effects are mean donor log2(CPM+1) differences, not raw-count fold changes','limits':'Cells are not replicates. Within-compartment fractions are not whole-tissue abundance. All provisional subclusters retained in broad-compartment DE; mixture remains a limitation. F3+F4 pooling adapts original disease contrast; not paired or longitudinal. No forced trajectory root.'}
(O/'analysis_design.json').write_text(json.dumps(manifest,indent=2));deps={k:importlib.util.find_spec(k) is not None for k in ['gseapy','cytotrace2','scanpy']};deps['Rscript']=shutil.which('Rscript')
if deps['Rscript']:
 r=subprocess.run([deps['Rscript'],'-e','cat("monocle=",requireNamespace("monocle",quietly=TRUE)," CytoTRACE=",requireNamespace("CytoTRACE",quietly=TRUE))'],capture_output=True,text=True,timeout=60);deps['R_packages']=r.stdout+r.stderr
(O/'trajectory_dependencies.json').write_text(json.dumps(deps,indent=2));results=[]
for task,labels in manual.items():
 try:
  record(task,status='running');g=task.split('_')[0];c=pd.read_csv(S/task/'subcluster_coordinates_metadata.tsv.gz',sep='\t',index_col=0,keep_default_na=False);c['broad_review_label']=c.subcluster_preliminary.astype(str).map(labels);c=c[c.disease_comparison_eligible.astype(str).str.lower().eq('true')].copy();c['group']=np.where(c.F_stage_from_Figure4C.astype(str).isin(['3','4','3.0','4.0']),'F3_F4','F0') if g=='GSE244832' else c['disease status'];c.to_csv(O/(task+'_eligible_cells.tsv.gz'),sep='\t');src=ad.read_h5ad(P/g/'qc_pass_counts.h5ad',backed='r');a=src[c.index,:].to_memory();src.file.close();a.obs=c
  for sorting,obs in c.groupby('sorting'):
   key=task+'_'+str(sorting).replace('+','positive').replace('-','negative');d=O/key;d.mkdir();tab=pd.crosstab(obs.subcluster_preliminary,obs.group);roe=tab/np.outer(tab.sum(axis=1),tab.sum(axis=0))*tab.to_numpy().sum();roe.to_csv(d/'subcluster_RoE.tsv',sep='\t');counts=pd.crosstab(obs.donor,obs.subcluster_preliminary);fractions=counts.div(counts.sum(axis=1),axis=0);fractions.join(obs.groupby('donor').group.first()).to_csv(d/'donor_subcluster_fractions.tsv',sep='\t');fig,ax=plt.subplots(figsize=(6,7));im=ax.imshow(roe,aspect='auto',vmin=0,vmax=2,cmap='RdBu_r');ax.set_xticks(range(len(roe.columns)),roe.columns);ax.set_yticks(range(len(roe)),[str(cl)+' '+labels[str(cl)] for cl in roe.index]);ax.set_title(key+' | descriptive Ro/e');fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(d/'F3CD_F4C_RoE.png',dpi=160);plt.close(fig)
   donors=sorted(obs.donor.unique());matrix=[];meta=[]
   for donor in donors:
    ci=np.flatnonzero((c.donor==donor)&(c.sorting==sorting));matrix.append(np.asarray(a.X[ci,:].sum(axis=0)).ravel());ds=c.iloc[ci];assert ds.group.nunique()==1;meta.append({'donor':donor,'group':ds.group.iloc[0],'cells':len(ci),'libraries':ds.library.nunique()})
   raw=np.vstack(matrix);meta=pd.DataFrame(meta);meta.to_csv(d/'independent_donor_manifest.tsv',sep='\t',index=False);groups=sorted(meta.group.unique());reference='F0' if g=='GSE244832' else next((z for z in groups if 'healthy' in z.lower()),None);case=next((z for z in groups if z!=reference),None)
   if not reference or not case:results.append({'task':key,'status':'comparison_unavailable','groups':groups});continue
   log=np.log2(raw/np.maximum(raw.sum(axis=1,keepdims=True),1)*1e6+1);xx=log[meta.group.eq(case)];yy=log[meta.group.eq(reference)];test=mannwhitneyu(xx,yy,axis=0,alternative='two-sided',method='asymptotic');effect=xx.mean(axis=0)-yy.mean(axis=0);symbols=a.var.gene_symbol_original.astype(str).to_numpy();de=pd.DataFrame({'gene_id':a.var_names,'gene_symbol':symbols,'mean_donor_log2CPM_difference':effect,'p_value':test.pvalue,'BH_q':bh(test.pvalue),'case_donors':len(xx),'reference_donors':len(yy)});de.to_csv(d/'all_gene_donor_rank_sum.tsv.gz',sep='\t',index=False);fig,ax=plt.subplots(figsize=(7,5));ax.scatter(effect,-np.log10(np.maximum(test.pvalue,1e-300)),s=3,alpha=.35);ax.axhline(-np.log10(.05),color='grey',ls='--');ax.set_xlabel('Mean donor log2(CPM+1) difference');ax.set_ylabel('-log10 unadjusted p');ax.set_title(key+' | '+case+' vs '+reference+' | donors '+str(len(xx))+'/'+str(len(yy)));fig.tight_layout();fig.savefig(d/'F3E_compartment_donor_DE.png',dpi=180);plt.close(fig);results.append({'task':key,'case':case,'reference':reference,'case_donors':len(xx),'reference_donors':len(yy),'genes':len(de),'raw_p_below_05':int((de.p_value<.05).sum()),'BH_q_below_05':int((de.BH_q<.05).sum()),'top_by_p':de.sort_values('p_value').head(5).to_dict('records')})
  del a;record(task,status='complete')
 except Exception:record(task,status='failed',error=traceback.format_exc())
(O/'results.json').write_text(json.dumps(results,indent=2));record('all',status='complete',comparisons=len(results))
