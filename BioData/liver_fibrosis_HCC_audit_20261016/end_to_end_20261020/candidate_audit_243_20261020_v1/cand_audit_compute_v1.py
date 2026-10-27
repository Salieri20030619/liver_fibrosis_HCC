from pathlib import Path
import json,hashlib,datetime,sys,tarfile,gzip,io,itertools,pickle
import numpy as np,pandas as pd
from scipy.stats import spearmanr,mannwhitneyu,ks_2samp
from scipy.spatial.distance import jensenshannon
from statsmodels.stats.multitest import multipletests
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC'
sys.path.insert(0,str(P/'python_deps'))
import joblib
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
O=P/'candidate_audit_243_20261020_v1';O.mkdir(exist_ok=False)
N=P/'continuation_newnetwork_models_v2';I=N/'features';C=P/'continuation_SCENIC_candidates_v1';S=P/'continuation_SCENIC_expanded_v1';M=P/'continuation_network_communication_inputs_v1/SCENIC'
def save(d,n):d.to_csv(O/n,sep='\t',index=False)
def js(x,n):(O/n).write_text(json.dumps(x,indent=2,default=str))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
X=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0);E=pd.read_csv(I/'external_features.tsv',sep='\t',index_col=0)
tm=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0).loc[X.index];em=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0).loc[E.index]
genes=X.columns.tolist();assert len(genes)==30 and genes==E.columns.tolist()
candidate=pd.read_csv(C/'F5G_candidates.tsv',sep='\t');assert set(candidate.TF)==set(genes)
selected=json.loads((N/'search/validation.json').read_text())['selected']['model'];assert selected=='Stack_RF_4__SVM_1.0'
inputs=[C/'F5G_candidates.tsv',I/'feature_provenance.json',I/'train_features.tsv',I/'external_features.tsv',I/'train_design.tsv',I/'external_design.tsv',I/'outer_cv_folds.tsv',I/'RF_4.joblib',I/'SVM_1.0.joblib',N/'search/meta_models.pkl',S/'regulons.csv',S/'adjacencies.tsv',M/'cell_metadata.tsv',P/'cont_SCENIC_postprocess.py',P/'cont_SCENIC_expanded.py',P/'cont_network_inputs_review_gated.py',P/'cont_newnetwork_generated_base.py',P/'cont_newnetwork_generated_combine.py',D/'HGNC_one_to_one_mapping_v1.tsv']
manifest=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in inputs]
save(pd.DataFrame(manifest),'frozen_input_manifest.tsv')
save(pd.DataFrame({'order':range(1,31),'gene':genes}),'frozen30_features.tsv')
js({'model':selected,'candidate_file_SHA256':sha(C/'F5G_candidates.tsv'),'genes':genes,'cohorts':[206,112],'endpoint':'F3/F4 versus F0/F1/F2 cross-sectional','main_workflow':'243 configurations, latest30 genes;113 optional supplement only','changes':'No model search, candidate replacement, endpoint change or external-label feature optimization'},'scope.json')
print('Frozen30',genes,flush=True)
# Model inspection and selector-only deterministic replay. No classifier retraining.
models={n:joblib.load(I/(n+'.joblib')) for n in ['RF_4','SVM_1.0']}
for mod in models.values():assert mod.n_features_in_==30
fold=pd.read_csv(I/'outer_cv_folds.tsv',sep='\t').set_index('sample').loc[X.index,'fold']
rferows=[]
for k in list(range(10))+['full']:
 use=np.ones(len(X),bool) if k=='full' else fold.to_numpy()!=k
 scale=StandardScaler().fit(X.to_numpy()[use]);z=scale.transform(X.to_numpy()[use])
 rfe=RFE(LogisticRegression(C=1,solver='liblinear',max_iter=2000),n_features_to_select=15,step=1).fit(z,tm.y.to_numpy()[use])
 for j,g in enumerate(genes):rferows.append(dict(fold=k,gene=g,selected=bool(rfe.support_[j]),rank=int(rfe.ranking_[j])))
 if k=='full':
  for mod in models.values():
   assert np.array_equal(rfe.support_,mod.named_steps['RFE'].support_)
   assert np.max(abs(scale.mean_-mod.named_steps['scale'].mean_))<1e-12
rfe=pd.DataFrame(rferows);save(rfe,'RFE_selector_only_replay.tsv')
with open(N/'search/meta_models.pkl','rb') as f:meta=pickle.load(f)[selected]
basepred=pd.DataFrame({n:mod.predict_proba(E.to_numpy())[:,1] for n,mod in models.items()},index=E.index)
pred=meta.predict_proba(basepred)[:,1];prior=pd.read_csv(N/'search/all243_external_predictions.tsv',sep='\t',index_col=0).loc[E.index,selected]
replay_error=float(np.max(abs(pred-prior)));assert replay_error<1e-10
js({'model_replay_max_error':replay_error,'full_RFE_saved_model_match':True,'fold_RFE':'reconstructed using exact original deterministic selector; historical fold selector objects were not persisted','full_selected_genes':list(X.columns[models['RF_4'].named_steps['RFE'].support_]),'RF_SVM_same_features':bool(np.array_equal(models['RF_4'].named_steps['RFE'].support_,models['SVM_1.0'].named_steps['RFE'].support_)),'meta_coefficients':meta.coef_.tolist(),'meta_intercept':meta.intercept_.tolist(),'RF_class_order':models['RF_4'].classes_.tolist(),'SVM_class_order':models['SVM_1.0'].classes_.tolist()},'model_replay.json')
jacc=[]
for k,l in itertools.combinations(range(10),2):
 a=set(rfe[(rfe.fold==k)&rfe.selected].gene);b=set(rfe[(rfe.fold==l)&rfe.selected].gene);jacc.append(dict(fold1=k,fold2=l,Jaccard=len(a&b)/len(a|b)))
save(pd.DataFrame(jacc),'RFE_fold_pair_stability.tsv')
print('RFE/model replay completed',flush=True)
# Raw counts, exact gene mapping and full normalization recomputation.
h=pd.read_csv(D/'HGNC_one_to_one_mapping_v1.tsv',sep='\t')
mapping=h.set_index('ensembl_gene_id').symbol;vectors={}
with tarfile.open(P/'GSE135251_RAW.tar') as tar:
 for mem in tar:
  gsm=mem.name.split('_')[0]
  if not mem.isfile() or gsm not in X.index:continue
  v=pd.read_csv(io.BytesIO(gzip.decompress(tar.extractfile(mem).read())),sep='\t',header=None,index_col=0).iloc[:,0]
  vectors[gsm]=v[~v.index.str.startswith('__')]
raws={'train':pd.DataFrame(vectors).loc[:,X.index],'external':pd.read_csv(D/'GSE162694_raw_counts.csv.gz',index_col=0).loc[:,E.index]}
mappingrows=[];counts={};librows=[];normerrors={}
for cohort,raw in raws.items():
 names=pd.Series(raw.index.str.split('.').str[0],index=raw.index).map(mapping)
 for g in genes:
  ids=names.index[names.eq(g)].tolist();mappingrows.append(dict(cohort=cohort,gene=g,raw_IDs=';'.join(ids),mapping_count=len(ids),HGNC_rows=int(h.symbol.eq(g).sum())))
  assert len(ids)==1
 counts[cohort]=pd.DataFrame({g:raw.loc[names.index[names.eq(g)][0]] for g in genes})
 totals=raw.sum(0);recomputed=np.log2(counts[cohort].div(totals,axis=0)*1e6+1)
 target=X if cohort=='train' else E;normerrors[cohort]=float(np.max(abs(recomputed.loc[target.index,genes].to_numpy()-target.to_numpy())));assert normerrors[cohort]<1e-10
 design=tm if cohort=='train' else em
 for sample,total in totals.items():librows.append(dict(cohort=cohort,sample=sample,F=int(design.loc[sample,'F']),y=int(design.loc[sample,'y']),library_count=float(total)))
 save(counts[cohort].rename_axis('sample').reset_index(),cohort+'_30_rawcounts.tsv')
save(pd.DataFrame(mappingrows),'raw_gene_mapping.tsv');save(pd.DataFrame(librows),'sample_library_counts.tsv')
js({'recomputed_log2CPM_max_abs_error':normerrors,'raw_shapes':{k:list(v.shape) for k,v in raws.items()},'raw_ID_duplicate_counts':{k:int(v.index.duplicated().sum()) for k,v in raws.items()},'train_external_exact_ID_intersection':list(set(X.index)&set(E.index)),'missing_feature_values':int(X.isna().sum().sum()+E.isna().sum().sum())},'mapping_validation.json')
# Descriptive effects; no feature elimination.
rows=[];frows=[]
for label,x,md in [('train',X,tm),('external',E,em)]:
 y=md.y.to_numpy();raw=counts[label]
 for g in genes:
  v=x[g].to_numpy();a=v[y==1];b=v[y==0];u,p=mannwhitneyu(a,b,alternative='two-sided');rho,sp=spearmanr(v,md.F)
  sd=np.sqrt(((len(a)-1)*np.var(a,ddof=1)+(len(b)-1)*np.var(b,ddof=1))/(len(a)+len(b)-2));d=(np.mean(a)-np.mean(b))/sd if sd>0 else np.nan
  rows.append(dict(cohort=label,gene=g,n=len(v),positive=int(y.sum()),detected_fraction=float((raw[g]>0).mean()),CPM_ge1_fraction=float((v>=1).mean()),mean=float(v.mean()),SD=float(v.std(ddof=1)),median=float(np.median(v)),q05=float(np.quantile(v,.05)),q95=float(np.quantile(v,.95)),mean_advanced=float(a.mean()),mean_early=float(b.mean()),mean_difference_advanced_minus_early=float(a.mean()-b.mean()),median_difference=float(np.median(a)-np.median(b)),cohens_d=float(d),rank_biserial=float(2*u/(len(a)*len(b))-1),unflipped_AUC=float(roc_auc_score(y,v)),MW_p=float(p),F_spearman=float(rho),F_spearman_p=float(sp)))
  for f,grp in md.groupby('F'):
   vv=x.loc[grp.index,g];frows.append(dict(cohort=label,gene=g,F=int(f),n=len(vv),mean=float(vv.mean()),median=float(vv.median()),SD=float(vv.std()),detected_fraction=float((raw.loc[grp.index,g]>0).mean())))
effects=pd.DataFrame(rows)
for label,ix in effects.groupby('cohort').groups.items():
 effects.loc[ix,'MW_BH_q30']=multipletests(effects.loc[ix,'MW_p'],method='fdr_bh')[1];effects.loc[ix,'F_spearman_BH_q30']=multipletests(effects.loc[ix,'F_spearman_p'],method='fdr_bh')[1]
save(effects,'gene_cohort_effects.tsv');save(pd.DataFrame(frows),'gene_per_F_distributions.tsv')
shifts=[]
for g in genes:
 for label,ti,ei in [('all',X.index,E.index)]+[('F'+str(f),tm.index[tm.F==f],em.index[em.F==f]) for f in range(5)]:
  if not len(ti) or not len(ei):continue
  a=X.loc[ti,g];b=E.loc[ei,g];ks,p=ks_2samp(a,b)
  shifts.append(dict(gene=g,stratum=label,train_n=len(a),external_n=len(b),external_minus_train_mean=float(b.mean()-a.mean()),shift_in_training_SD=float((b.mean()-a.mean())/X[g].std()),KS=float(ks),KS_p=float(p),external_below_train_min=float((b<a.min()).mean()),external_above_train_max=float((b>a.max()).mean())))
sh=pd.DataFrame(shifts);save(sh,'cross_cohort_feature_shifts.tsv')
corr=[]
for cohort,x in [('train',X),('external',E)]:
 c=x.corr(method='spearman')
 for g,j in itertools.combinations(genes,2):corr.append(dict(cohort=cohort,gene1=g,gene2=j,rho=c.loc[g,j]))
save(pd.DataFrame(corr),'feature_pair_correlations.tsv')
# Upstream cell/sample composition and stratified RSS diagnostics.
m=pd.read_csv(M/'cell_metadata.tsv',sep='\t',index_col=0);a=pd.read_csv(S/'aucell_matrix.csv',index_col=0);assert a.index.equals(m.index)
allmeta=pd.read_csv(P/'continuation_annotation_refined_v1/cell_metadata_reviewed.tsv.gz',sep='\t',index_col=0)
my=allmeta[allmeta.reviewed_broad_lineage.eq('Myeloid')].copy();my['source_record']=my.study+'|'+my.donor.astype(str)
record=m.groupby('source_record').agg(study=('study','first'),author_donor=('donor','first'),GSMs=('GSM',lambda z:';'.join(sorted(set(z)))),F=('F','first'),F_unique=('F','nunique'),sampled_cells=('F','size'),states=('refined_cluster','nunique'))
record['all_myeloid_cells']=my.groupby('source_record').size();record['independent_donor_status']='source-scoped label; cross-study identity unconfirmed';save(record.reset_index(),'network54_records.tsv')
comp=m.groupby(['study','source_record','F','refined_cluster','regulon_role']).size().rename('sampled_cells').reset_index()
full=my.groupby(['study','source_record','F','refined_cluster']).size().rename('all_myeloid_cells').reset_index();comp=comp.merge(full,how='left');save(comp,'sampled_record_state_composition.tsv')
save(m.groupby(['study','F','refined_cluster']).size().rename('cells').reset_index(),'network_study_stage_state_cells.tsv')
save(record.reset_index().groupby(['study','F']).size().rename('source_records').reset_index(),'network_study_stage_records.tsv')
dis=m.regulon_role.eq('discovery_advanced_F3_F4');md=m.loc[dis];ad=a.loc[dis];good=ad.std().gt(1e-14);ad=ad.loc[:,good];states=sorted(md.refined_cluster.unique())
def rss(av,labels,weights=None):
 weights=np.ones(len(av)) if weights is None else np.asarray(weights)
 arr=av.to_numpy()*weights[:,None];arr=arr/arr.sum(0)
 rr=[]
 for state in sorted(set(labels)):
  t=(np.asarray(labels)==state)*weights;t=t/t.sum()
  vals=1-jensenshannon(arr,t[:,None],axis=0)
  ranks=pd.Series(vals,index=av.columns).rank(ascending=False,method='first')
  rr.extend(dict(state=state,regulon=g,RSS=float(v),rank=int(ranks[g])) for g,v in zip(av.columns,vals))
 return pd.DataFrame(rr)
orig=rss(ad,md.refined_cluster)
orig_saved=pd.read_csv(C/'all_RSS.tsv',sep='\t');check=orig.merge(orig_saved,left_on=['state','regulon'],right_on=['subcluster','regulon']);assert np.max(abs(check.RSS_x-check.RSS_y))<1e-10
rsn=md.groupby(['source_record','refined_cluster']).F.transform('size');rn=md.groupby('source_record').F.transform('size')
variants={'original_cells':orig,'record_state_weighted_cells':rss(ad,md.refined_cluster,1/rsn),'record_weighted_cells':rss(ad,md.refined_cluster,1/rn)}
agg=ad.groupby([md.source_record,md.refined_cluster]).mean();variants['record_state_aggregated']=rss(agg,agg.index.get_level_values(1))
for study in sorted(md.study.unique()):
 ix=md.study.eq(study);variants['study_only_'+study]=rss(ad.loc[ix],md.loc[ix,'refined_cluster'])
allrss=[];overlap=[]
for name,r in variants.items():
 r=r.assign(variant=name);allrss.append(r)
 for state,g in r.groupby('state'):
  old=set(orig[(orig.state==state)&(orig['rank']<=10)].regulon);new=set(g[g['rank']<=10].regulon)
  overlap.append(dict(variant=name,state=state,overlap_top10=len(old&new),Jaccard=len(old&new)/len(old|new)))
save(pd.concat(allrss),'RSS_weighting_diagnostics.tsv');save(pd.DataFrame(overlap),'RSS_top10_overlap_decomposition.tsv')
# Leave one discovery source record out; ranking audit only, no new candidates.
loo=[]
for recordid in md.source_record.unique():
 ix=md.source_record.ne(recordid);r=rss(ad.loc[ix],md.loc[ix,'refined_cluster'])
 for _,c in candidate.iterrows():
  row=r[(r.state==c.state)&(r.regulon==c.regulon)]
  loo.append(dict(left_out_record=recordid,gene=c.TF,state=c.state,rank=float(row['rank'].iloc[0]) if len(row) else np.nan,top10=bool(len(row) and row['rank'].iloc[0]<=10)))
save(pd.DataFrame(loo),'candidate_RSS_leave_one_record_out.tsv')
edges=pd.read_csv(C/'all_regulon_targets.tsv',sep='\t');edges[edges.TF.isin(genes)].to_csv(O/'candidate_regulon_targets.tsv',sep='\t',index=False)
print('Bulk, mapping, RSS completed',flush=True)
# One row per candidate, preserving evidence layers.
imp=pd.read_csv(N/'SHAP/gene_importance.tsv',sep='\t',index_col=0).iloc[:,0]
merged=[]
for g in genes:
 c=candidate[candidate.TF==g].iloc[0];tr=effects[(effects.cohort=='train')&(effects.gene==g)].iloc[0];ex=effects[(effects.cohort=='external')&(effects.gene==g)].iloc[0]
 rf=rfe[(rfe.gene==g)&(rfe.fold!='full')];ll=pd.DataFrame(loo);ll=ll[ll.gene==g]
 row=dict(gene=g,module=int(c.module),discovery_state=c.state,regulon=c.regulon,original_RSS=float(c.RSS),original_RSS_rank=int(c['rank']),selection_rule='module-highest AUCell state intersect state RSS top10; not fibrosis-DE filter',regulon_target_count=int(edges.TF.eq(g).sum()),train_mapping_n=int(next(z['mapping_count'] for z in mappingrows if z['cohort']=='train' and z['gene']==g)),external_mapping_n=int(next(z['mapping_count'] for z in mappingrows if z['cohort']=='external' and z['gene']==g)),RFE_outer_selected_count=int(rf.selected.sum()),RFE_outer_denominator=10,RFE_full_selected=bool(rfe[(rfe.gene==g)&(rfe.fold=='full')].selected.iloc[0]),mean_abs_SHAP=float(imp[g]),SHAP_rank=int(imp.rank(ascending=False,method='min')[g]),RSS_leave_one_record_top10_fraction=float(ll.top10.mean()),RSS_leave_one_record_worst_rank=float(ll['rank'].max()),direction_consistent=bool(np.sign(tr.mean_difference_advanced_minus_early)==np.sign(ex.mean_difference_advanced_minus_early)),external_shift_training_SD=float(sh[(sh.gene==g)&(sh.stratum=='all')].shift_in_training_SD.iloc[0]))
 for prefix,er in [('train',tr),('external',ex)]:
  for key in ['detected_fraction','CPM_ge1_fraction','SD','mean_difference_advanced_minus_early','cohens_d','rank_biserial','unflipped_AUC','MW_p','MW_BH_q30','F_spearman','F_spearman_BH_q30']:row[prefix+'_'+key]=er[key]
 for name in variants:
  rr=variants[name];z=rr[(rr.state==c.state)&(rr.regulon==c.regulon)];row['RSS_rank_'+name]=int(z['rank'].iloc[0]) if len(z) else np.nan
 row['evidence_level']='coexpression+motif-supported computational regulon; descriptive bulk association; conditional model contribution; no direct binding or causal validation'
 merged.append(row)
save(pd.DataFrame(merged),'30_feature_audit.tsv')
summary={'status':'completed','scope':'candidate audit only; mainmodel unchanged;113 optional supplementary','genes':genes,'source_records':len(record),'discovery_records':int(md.source_record.nunique()),'sampled_cells':len(m),'discovery_cells':len(md),'discovery_states':states,'RFE_pairwise_Jaccard_median':float(pd.DataFrame(jacc).Jaccard.median()),'train_BH05_genes':effects[(effects.cohort=='train')&(effects.MW_BH_q30<.05)].gene.tolist(),'external_BH05_genes':effects[(effects.cohort=='external')&(effects.MW_BH_q30<.05)].gene.tolist(),'same_direction_genes':int(pd.DataFrame(merged).direction_consistent.sum()),'mapping_max_error':normerrors,'model_replay_error':replay_error,'original_input_hashes_unchanged':all(sha(Path(r['path']))==r['sha256'] for r in manifest)}
js(summary,'validation.json');print(json.dumps(summary),flush=True)
