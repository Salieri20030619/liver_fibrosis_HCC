from pathlib import Path
import sys,json,pickle,hashlib,time,datetime
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));O=P/'independent_validation_A_20261021_v1';N=P/'continuation_newnetwork_models_v2';I=N/'features';S=N/'search';D=P.parent/'analysis_stage1_20261017_075157UTC'
import numpy as np,pandas as pd,joblib
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss,confusion_matrix,roc_curve
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def js(x,n):(O/n).write_text(json.dumps(x,indent=2,default=str))
assert not (O/'prediction_metrics.tsv').exists(),'Do not rerun completed validation'
t0=time.time();X=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0);tm=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0).loc[X.index];oldE=pd.read_csv(I/'external_features.tsv',sep='\t',index_col=0);em=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0).loc[oldE.index];models={n:joblib.load(I/(n+'.joblib')) for n in ['RF_4','SVM_1.0']}
with open(S/'meta_models.pkl','rb') as f:meta=pickle.load(f)['Stack_RF_4__SVM_1.0']
mask=models['RF_4'].named_steps['RFE'].support_;assert np.array_equal(mask,models['SVM_1.0'].named_steps['RFE'].support_);active=list(X.columns[mask]);files=[I/'RF_4.joblib',I/'SVM_1.0.joblib',S/'meta_models.pkl',I/'train_features.tsv',I/'external_features.tsv',S/'all243_training_OOF.tsv',S/'all243_external_predictions.tsv',N/'validation/fixed_threshold_performance.tsv'];manifest={str(p):sha(p) for p in files};js(manifest,'frozen_model_hashes.json')
raw=pd.read_csv(O/'GSE174478_counts.txt.gz',sep='\t',index_col=0);assert raw.columns.is_unique and raw.shape[1]==94;assert np.isfinite(raw).all().all() and raw.to_numpy().min()>=0 and (raw.sum()>0).all()
dup=raw.loc[raw.index.duplicated(False)];dupcheck=[]
for gene,z in dup.groupby(level=0):dupcheck.append(dict(ensembl=gene,rows=len(z),identical_values=bool((z.to_numpy()==z.iloc[0].to_numpy()).all())))
pd.DataFrame(dupcheck).to_csv(O/'duplicate_gene_audit.tsv',sep='\t',index=False);dup_mass=raw.loc[raw.index.duplicated(keep='first')].sum()/raw.sum();js({'duplicate_IDs':len(dupcheck),'duplicate_rows':len(dup),'extra_copy_mass_fraction_min':float(dup_mass.min()),'extra_copy_mass_fraction_max':float(dup_mass.max()),'handling':'Retain original full deposited column sum as frozen plan; exclude ambiguous mapped symbols from feature matrix; no deduplication or altered CPM denominator','first_attempt':'Overstrict raw unique-ID gate stopped beforeprediction; preserved worker and log; this version uses already declared ambiguous-mapping exclusion'},'duplicate_handling.json')
m=pd.read_csv(O/'GSE174478_sample_metadata.tsv',sep='\t');m['sample']=m.title.str.replace('Liver transcriptome of ','',regex=False);assert m['sample'].is_unique and m.GSM.is_unique and set(m['sample'])==set(raw.columns);m=m.set_index('sample').loc[raw.columns];m['F']=pd.to_numeric(m['fibrosis stage']);assert m.F.isin(range(5)).all();m['y']=(m.F>=3).astype(int);assert (m.tissue=='liver').all() and (m['disease state']=='NAFLD').all();m[['GSM','title','F','y','age','Sex','instrument_model']].to_csv(O/'new_sample_alignment.tsv',sep='\t')
# Check known public accession overlap, not genetic donor deduplication.
new_gsms=set(m.GSM);known=set()
for frame in [tm,em]:
 known.update(v for v in frame.index.astype(str) if v.startswith('GSM'))
 for col in frame:
  known.update(v for v in frame[col].astype(str) if v.startswith('GSM'))
assert not new_gsms&known
h=pd.read_csv(D/'HGNC_one_to_one_mapping_v1.tsv',sep='\t').set_index('ensembl_gene_id').symbol;ids=raw.index.astype(str);mapped=pd.Series(ids.str.split('.').str[0],index=ids).map(h);valid=mapped.notna()&~mapped.duplicated(keep=False);norm=np.log2(raw.div(raw.sum(0),axis=1)*1e6+1).loc[valid.to_numpy()];norm.index=mapped[valid].values;assert norm.index.is_unique
missing=sorted(set(X.columns)-set(norm.index));assert not set(active)&set(missing),'Active15 incomplete; stop before prediction'
E=norm.reindex(X.columns).T
for g in missing:E[g]=float(X[g].mean())
assert E.index.equals(m.index) and list(E.columns)==list(X.columns) and np.isfinite(E).all().all();E.to_csv(O/'new_features_original_scale.tsv',sep='\t')
pd.DataFrame({'gene':X.columns,'active15':mask,'available':X.columns.isin(norm.index),'imputed_training_mean':X.columns.isin(missing)}).to_csv(O/'feature_coverage.tsv',sep='\t',index=False)
qual={'utc_before_predictions':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selected':'GSE174478','samples':len(m),'cases':int(m.y.sum()),'controls':int((m.y==0).sum()),'F_counts':m.F.value_counts().sort_index().to_dict(),'raw_dimensions':list(raw.shape),'integer_values':bool(np.equal(raw,np.floor(raw)).all().all()),'mapped_genes':len(norm),'active_coverage':15,'input30_coverage':30-len(missing),'inactive_missing_imputed':missing,'GSM_overlap_known':list(new_gsms&known),'independence':'Tokyo Japanese collection versus European train and US old external; no known duplicateGSM; not genotype confirmed','preprocessing':'Original full deposited gene-count column sum, log2(CPM+1), fixed HGNCone-to-one; training scaler/RFESupport unchanged; no recalibration','counts_sha256':sha(O/'GSE174478_counts.txt.gz'),'metadata_sha256':sha(O/'GSE174478_sample_metadata.tsv')};js(qual,'QUALIFIED_BEFORE_PREDICTION.json')
print(json.dumps(qual),flush=True)
def predict(z):
 a=np.column_stack([models[n].predict_proba(z.to_numpy())[:,1] for n in ['RF_4','SVM_1.0']]);return meta.predict_proba(a)[:,1]
oldpred=predict(oldE);saved=pd.read_csv(S/'all243_external_predictions.tsv',sep='\t',index_col=0).loc[oldE.index,'Stack_RF_4__SVM_1.0'];err=float(np.max(abs(oldpred-saved)));assert err<1e-12
pp=predict(E);pd.DataFrame({'GSM':m.GSM,'F':m.F,'y':m.y,'probability':pp},index=E.index).to_csv(O/'new_predictions.tsv',sep='\t')
threshold=float(pd.read_csv(N/'validation/fixed_threshold_performance.tsv',sep='\t').fixed_training_Youden_threshold.iloc[0]);oof=pd.read_csv(S/'all243_training_OOF.tsv',sep='\t',index_col=0).loc[X.index,'Stack_RF_4__SVM_1.0'];fpr,tpr,thr=roc_curve(tm.y,oof);assert abs(threshold-thr[np.argmax(tpr-fpr)])<1e-12
rows=[];shifts=[];outliers=[]
for cohort,z,yy,pr,design in [('GSE162694',oldE,em.y.to_numpy(),oldpred,em),('GSE174478',E,m.y.to_numpy(),pp,m)]:
 rng=np.random.default_rng(20261021);pos=np.flatnonzero(yy==1);neg=np.flatnonzero(yy==0);indices=np.array([np.r_[rng.choice(pos,len(pos)),rng.choice(neg,len(neg))] for j in range(2000)]);np.save(O/(cohort+'_bootstrap_indices.npy'),indices);vals=[roc_auc_score(yy[ix],pr[ix]) for ix in indices];lo,hi=np.quantile(vals,[.025,.975]);tn,fp,fn,tp=confusion_matrix(yy,pr>=threshold).ravel();rows.append(dict(cohort=cohort,n=len(yy),cases=len(pos),controls=len(neg),AUC=roc_auc_score(yy,pr),CI_low=lo,CI_high=hi,includes_0_5=bool(lo<=.5<=hi),AP=average_precision_score(yy,pr),Brier=brier_score_loss(yy,pr),threshold=threshold,TN=int(tn),FP=int(fp),FN=int(fn),TP=int(tp),sensitivity=tp/(tp+fn),specificity=tn/(tn+fp),PPV=tp/(tp+fp) if tp+fp else None,NPV=tn/(tn+fn) if tn+fn else None))
 zz=models['RF_4'].named_steps['scale'].transform(z.to_numpy())[:,mask];hits=abs(zz)>3;outliers.append(dict(cohort=cohort,features=15,numerator=int(hits.sum()),denominator=int(hits.size),fraction=float(hits.mean()),samples_any=int(hits.any(1).sum())))
 for fs in ['all',0,1,2,3,4]:
  tx=X if fs=='all' else X.loc[tm.F==fs];ex=z if fs=='all' else z.loc[design.F==fs]
  for g in active:shifts.append(dict(cohort=cohort,F=fs,gene=g,train_n=len(tx),external_n=len(ex),external_minus_training_mean=float(ex[g].mean()-tx[g].mean()),difference_in_training_SD=float((ex[g].mean()-tx[g].mean())/X[g].std(ddof=0))))
pd.DataFrame(rows).to_csv(O/'prediction_metrics.tsv',sep='\t',index=False);pd.DataFrame(shifts).to_csv(O/'expression_shift_by_stage.tsv',sep='\t',index=False);pd.DataFrame(outliers).to_csv(O/'active15_outliers.tsv',sep='\t',index=False)
assert all(sha(Path(p))==h for p,h in manifest.items());js({'status':'completed','model_fit_calls':0,'prediction_replay_old_max_error':err,'new_replay_max_error':float(np.max(abs(predict(E)-pp))),'hashes_unchanged':True,'selected_cohorts_predicted':1,'threshold':threshold,'seconds':time.time()-t0,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},'validation.json');print(pd.DataFrame(rows).to_string(index=False),flush=True)
