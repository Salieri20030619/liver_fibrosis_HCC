from pathlib import Path
import sys,json,hashlib,time,tarfile,gzip,io,pickle
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'A_cohort_shift_diagnostic_20261021_v1';O.mkdir(exist_ok=False);I=P/'continuation_newnetwork_models_v2/features';S=P/'continuation_newnetwork_models_v2/search'
import pandas as pd,numpy as np,joblib
def js(x,n):(O/n).write_text(json.dumps(x,indent=2,default=str))
def tab(x,n):pd.DataFrame(x).to_csv(O/n,sep='\t',index=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
start=time.time();X=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0);E=pd.read_csv(I/'external_features.tsv',sep='\t',index_col=0);tm=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0).loc[X.index];em=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0).loc[E.index];genes=X.columns.tolist();assert E.columns.tolist()==genes
sources=[I/'train_features.tsv',I/'external_features.tsv',I/'train_design.tsv',I/'external_design.tsv',I/'RF_4.joblib',I/'SVM_1.0.joblib',S/'meta_models.pkl',S/'all243_external_predictions.tsv'];manifest={str(p):sha(p) for p in sources};js(manifest,'baseline_source_hashes.json');js({'scope':'A only read-only diagnosis, no training or correction','scale':'existing log2(count/full gene-assigned sum*1e6+1)','uncertainty':'2000 within-F resamples seed20261021; pooled early/late weighted by observed stage counts; standardized effects use equal weights F0/F1/F2 and F3/F4 separately; same-F cohort differences; percentile CIs descriptive, no multiplicity claim','global_background':'all common version-stripped Ensembl IDs and separately genes CPM>=1 in at least10% of samples in BOTH cohorts; diagnostic only, not feature selection','limits':'No external-label calibration, correction, feature selection or causal identification'},'frozen_diagnostic_plan.json')
models={n:joblib.load(I/(n+'.joblib')) for n in ['RF_4','SVM_1.0']};mask=models['RF_4'].named_steps['RFE'].support_;assert np.array_equal(mask,models['SVM_1.0'].named_steps['RFE'].support_)
for m in models.values():assert np.max(abs(m.named_steps['scale'].mean_-X.mean().to_numpy()))<1e-12;assert np.max(abs(m.named_steps['scale'].scale_-X.std(ddof=0).to_numpy()))<1e-12
with open(S/'meta_models.pkl','rb') as f:meta=pickle.load(f)['Stack_RF_4__SVM_1.0']
pred=meta.predict_proba(pd.DataFrame({n:m.predict_proba(E.to_numpy())[:,1] for n,m in models.items()},index=E.index))[:,1];old=pd.read_csv(S/'all243_external_predictions.tsv',sep='\t',index_col=0).loc[E.index,'Stack_RF_4__SVM_1.0'];prederr=float(np.max(abs(pred-old)))
assert X.index.is_unique and E.index.is_unique and np.isfinite(X).all().all() and np.isfinite(E).all().all()
for md in [tm,em]:assert np.array_equal(md.y.to_numpy(),(md.F>=3).astype(int).to_numpy())
vectors={};special=[]
with tarfile.open(P/'GSE135251_RAW.tar') as tar:
 for mem in tar:
  gsm=mem.name.split('_')[0]
  if not mem.isfile() or gsm not in X.index:continue
  v=pd.read_csv(io.BytesIO(gzip.decompress(tar.extractfile(mem).read())),sep='\t',header=None,index_col=0).iloc[:,0];special.append(dict(sample=gsm,excluded_HTSeq_summary_count=float(v[v.index.str.startswith('__')].sum())));vectors[gsm]=v[~v.index.str.startswith('__')]
raw={'train':pd.DataFrame(vectors).loc[:,X.index],'external':pd.read_csv(D/'GSE162694_raw_counts.csv.gz',index_col=0).loc[:,E.index]};tab(special,'HTSeq_summary_excluded.tsv');h=pd.read_csv(D/'HGNC_one_to_one_mapping_v1.tsv',sep='\t');mapping=h.set_index('ensembl_gene_id').symbol;norm={};units=[];lib=[];maprows=[];low=[];errors={}
for label,r in raw.items():
 assert r.index.is_unique and np.isfinite(r).all().all() and (r>=0).all().all();r.index=r.index.str.split('.').str[0];assert r.index.is_unique;tot=r.sum();norm[label]=np.log2(r.div(tot,axis=1)*1e6+1);target=X if label=='train' else E;md=tm if label=='train' else em;names=pd.Series(r.index,index=r.index).map(mapping)
 ids=[]
 for g in genes:
  gi=names.index[names.eq(g)].tolist();assert len(gi)==1;ids+=gi;maprows.append(dict(cohort=label,gene=g,Ensembl=gi[0],mapping_n=1))
  counts=r.loc[gi[0]];low.append(dict(cohort=label,gene=g,zero_fraction=float((counts==0).mean()),CPM_below1_fraction=float((counts/tot*1e6<1).mean()),missing_n=int(counts.isna().sum())))
 z=norm[label].loc[ids].T;z.columns=genes;errors[label]=float(np.max(abs(z.to_numpy()-target.to_numpy())));assert errors[label]<1e-10
 units.append(dict(cohort=label,rows=r.shape[0],samples=r.shape[1],noninteger_fraction=float((abs(r.to_numpy()-np.rint(r.to_numpy()))>1e-6).mean()),max_value=float(r.max().max()),zero_fraction=float((r==0).to_numpy().mean()),library_min=float(tot.min()),library_median=float(tot.median()),library_max=float(tot.max()),mapped_rows=int(names.notna().sum())))
 for sample in r:
  lib.append(dict(cohort=label,sample=sample,F=int(md.loc[sample,'F']),total_gene_counts=float(tot[sample]),top10_gene_count_fraction=float(r[sample].nlargest(10).sum()/tot[sample])))
tab(units,'raw_input_numeric_semantics.tsv');tab(lib,'library_composition.tsv');tab(maprows,'candidate_mapping.tsv');tab(low,'candidate_zero_low_expression.tsv');print('Raw semantics and normalization checked',flush=True)
common=raw['train'].index.intersection(raw['external'].index);background=[]
bgshift=norm['external'].loc[common].mean(axis=1)-norm['train'].loc[common].mean(axis=1);expressed=((norm['train'].loc[common]>=1).mean(axis=1)>=.1)&((norm['external'].loc[common]>=1).mean(axis=1)>=.1)
for label,r in raw.items():
 frac=r.loc[common].sum()/r.sum();tab([dict(sample=s,common_gene_count_fraction=float(v)) for s,v in frac.items()],label+'_common_denominator_mass.tsv')
for label,idx in [('all_common',common),('CPM1_at_least10pct_both',common[expressed])]:
 for f in ['pooled',0,1,2,3,4]:
  a=norm['train'].loc[idx] if f=='pooled' else norm['train'].loc[idx,tm.index[tm.F==f]];b=norm['external'].loc[idx] if f=='pooled' else norm['external'].loc[idx,em.index[em.F==f]];v=b.mean(axis=1)-a.mean(axis=1);background.append(dict(background=label,F=f,genes=len(idx),train_n=a.shape[1],external_n=b.shape[1],median_shift=float(v.median()),q10=float(v.quantile(.1)),q90=float(v.quantile(.9)),positive_fraction=float((v>0).mean())))
tab(background,'global_background_shift.tsv');tab(pd.DataFrame({'Ensembl':common,'gene':pd.Series(common).map(mapping),'external_minus_train_mean':bgshift.to_numpy(),'expressed_both':expressed.to_numpy()}),'global_gene_shifts.tsv')
# Bootstrap stratified by exact F; no modification of predictors.
rng=np.random.default_rng(20261021);means={};boots={};counts={};group=[]
for label,x,md in [('train',X,tm),('external',E,em)]:
 means[label]={};boots[label]={};counts[label]={}
 for f in range(5):
  a=x.loc[md.F==f].to_numpy();counts[label][f]=len(a);means[label][f]=a.mean(axis=0);ix=rng.integers(0,len(a),(2000,len(a)));boots[label][f]=a[ix].mean(axis=1)
  for j,g in enumerate(genes):
   lo,hi=np.quantile(boots[label][f][:,j],[.025,.975]);group.append(dict(cohort=label,F=f,n=len(a),gene=g,mean=means[label][f][j],CI_low=lo,CI_high=hi))
tab(group,'same_F_group_means_CI.tsv');effects={};effectboot={};rows=[]
for label in means:
 for method in ['observed_composition','equal_stage_composition']:
  weights={f:(counts[label][f]/sum(counts[label][k] for k in ([0,1,2] if f<3 else [3,4])) if method=='observed_composition' else (1/3 if f<3 else 1/2)) for f in range(5)}
  delta=sum((1 if f>=3 else -1)*weights[f]*means[label][f] for f in range(5));bb=sum((1 if f>=3 else -1)*weights[f]*boots[label][f] for f in range(5));effects[label,method]=delta;effectboot[label,method]=bb
  for j,g in enumerate(genes):
   lo,hi=np.quantile(bb[:,j],[.025,.975]);rows.append(dict(cohort=label,method=method,gene=g,advanced_minus_early=delta[j],CI_low=lo,CI_high=hi,includes_zero=bool(lo<=0<=hi)))
tab(rows,'within_cohort_effects_CI.tsv');rev=[];shifts=[]
for j,g in enumerate(genes):
 row={'gene':g,'used15':bool(mask[j])}
 for method in ['observed_composition','equal_stage_composition']:
  a=effects['train',method][j];b=effects['external',method][j];la,ha=np.quantile(effectboot['train',method][:,j],[.025,.975]);lb,hb=np.quantile(effectboot['external',method][:,j],[.025,.975]);row.update({method+'_train_delta':a,method+'_external_delta':b,method+'_sign_reversal':bool(a*b<0),method+'_both_CI_opposite_nonzero':bool((ha<0 and lb>0)or(hb<0 and la>0))})
 rev.append(row)
 for f in range(5):
  delta=means['external'][f][j]-means['train'][f][j];bb=boots['external'][f][:,j]-boots['train'][f][:,j];lo,hi=np.quantile(bb,[.025,.975]);shifts.append(dict(gene=g,F=f,train_n=counts['train'][f],external_n=counts['external'][f],external_minus_train=delta,in_training_SD=delta/X[g].std(ddof=0),CI_low=lo,CI_high=hi))
tab(rev,'direction_reversal_diagnostics.tsv');tab(shifts,'same_F_cross_cohort_shift_CI.tsv');zr=[];zsum={}
for label,x in [('train',X),('external',E)]:
 z=models['RF_4'].named_steps['scale'].transform(x.to_numpy());hit=abs(z)>3;used=hit[:,mask];zsum[label]={'numerator':int(used.sum()),'denominator':int(used.size),'fraction':float(used.mean()),'samples_with_any':int(used.any(axis=1).sum()),'samples':len(x),'features':int(mask.sum())}
 for j,g in enumerate(genes):zr.append(dict(cohort=label,gene=g,used15=bool(mask[j]),outlier_count=int(hit[:,j].sum()),sample_denominator=len(x),outlier_fraction=float(hit[:,j].mean()),mean_z=float(z[:,j].mean())))
tab(zr,'frozen_scaler_outlier_by_gene.tsv');js(zsum,'frozen_scaler_outlier_denominators.json');assert all(sha(Path(p))==h for p,h in manifest.items());v={'status':'completed','seconds':time.time()-start,'normalization_replay_error':errors,'prediction_replay_error':prederr,'no_fit_calls':True,'source_hashes_unchanged':True,'F_counts':counts,'raw_gene_rows':{k:len(r) for k,r in raw.items()},'common_Ensembl':len(common),'direction_reversals':sum(r['observed_composition_sign_reversal'] for r in rev),'opposite_nonzero_CIs':sum(r['observed_composition_both_CI_opposite_nonzero'] for r in rev),'standardized_reversals':sum(r['equal_stage_composition_sign_reversal'] for r in rev),'outliers':zsum};js(v,'validation.json');print(json.dumps(v),flush=True)
