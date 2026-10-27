from pathlib import Path
import sys,json,tarfile,gzip,io,hashlib,traceback,itertools
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'))
import numpy as np,pandas as pd
from sklearn.model_selection import StratifiedKFold,cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier,StackingClassifier
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import roc_auc_score,brier_score_loss
import joblib
D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'F8_F9_real_full_cohort_v1';O.mkdir(exist_ok=False)
try:
 rows=[];vectors=[];genes=None
 clinical={z['GSM']:z for z in json.loads((P/'GSE135251_parsed_samples.json').read_text())}
 with tarfile.open(P/'GSE135251_RAW.tar') as tar:
  for member in tar:
   if not member.isfile():continue
   gsm=member.name.split('_')[0];meta=clinical[gsm];c=dict(v.split(': ',1) for v in meta['characteristics_ch1'])
   if c['disease']!='NAFLD':continue
   f=int(c['fibrosis stage']);v=pd.read_csv(io.BytesIO(gzip.decompress(tar.extractfile(member).read())),sep='\t',header=None,index_col=0).iloc[:,0]
   v=v[~v.index.str.startswith('__')];assert v.index.is_unique
   if genes is None:genes=v.index
   assert genes.equals(v.index)
   vectors.append(v.to_numpy());rows.append({'sample':gsm,'donor':meta['title'][0],'F':f,'y':int(f>=3),'cohort':'GSE135251'})
 train=pd.DataFrame(vectors,index=[v['sample'] for v in rows],columns=genes).T;tm=pd.DataFrame(rows).set_index('sample');assert len(tm)==206 and tm.donor.is_unique
 ext=pd.read_csv(D/'GSE162694_raw_counts.csv.gz',index_col=0);em=pd.read_csv(D/'expression_metadata_alignment.tsv',sep='\t',keep_default_na=False).set_index('expression_id').loc[ext.columns]
 em=em[em['fibrosis stage'].isin(['0','1','2','3','4'])].copy();ext=ext[em.index];em['F']=em['fibrosis stage'].astype(int);em['y']=(em.F>=3).astype(int);assert len(em)==112
 h=pd.read_csv(D/'HGNC_one_to_one_mapping_v1.tsv',sep='\t').set_index('ensembl_gene_id').symbol
 candidates=pd.read_csv(P/'F5_to_F6_original_direction_v2/F5G_candidates.tsv',sep='\t');tf=set(candidates.TF)
 def expression(x):
  ids=x.index.str.split('.').str[0];assert ids.is_unique
  x.index=ids;norm=np.log2(x.div(x.sum(axis=0),axis=1)*1e6+1);names=pd.Series(ids,index=ids).map(h);norm=norm.loc[names.notna()];norm.index=names.dropna();assert norm.index.is_unique;return norm
 tr=expression(train);ex=expression(ext);features=sorted(tf&set(tr.index)&set(ex.index));features=[f for f in features if tr.loc[f].std()>0];assert len(features)>=3
 X=tr.loc[features].T.to_numpy();E=ex.loc[features].T.to_numpy();y=tm.y.to_numpy();ey=em.y.to_numpy()
 tm.to_csv(O/'train_design.tsv',sep='\t');em.to_csv(O/'external_design.tsv',sep='\t');tr.loc[features].T.to_csv(O/'train_features.tsv',sep='\t');ex.loc[features].T.to_csv(O/'external_features.tsv',sep='\t')
 (O/'feature_provenance.json').write_text(json.dumps({'genes':features,'source':'F5 discovery-only module/RSS intersection; no selection on bulk outcomes','candidate_sha256':hashlib.sha256((P/'F5_to_F6_original_direction_v2/F5G_candidates.tsv').read_bytes()).hexdigest(),'endpoint':'cross-sectional F3-4 versus F0-2; not incident HCC','original_gap':'Paper initial feature pool and exact113 grids undocumented; this is a prespecified reconstructed candidate interface, not exact replication','identity':'GSE135251 donor titles unique; GSE162694 expression IDs unique, no separate patient identifier. No known source reuse; cross-study genetic independence unverified.'},indent=2))
 print('PREPARED',X.shape,E.shape,'F counts',tm.F.value_counts().to_dict(),em.F.value_counts().to_dict(),flush=True)
 cv=StratifiedKFold(10,shuffle=True,random_state=20261020);folds=np.zeros(len(y),int)
 for k,(_,test) in enumerate(cv.split(X,y)):folds[test]=k
 pd.DataFrame({'sample':tm.index,'fold':folds}).to_csv(O/'outer_cv_folds.tsv',sep='\t',index=False)
 estimators={}
 for c in [.1,1.0]:
  for name,penalty in [('Lasso','l1'),('Ridge','l2'),('ElasticNet','elasticnet')]:estimators[name+'_'+str(c)]=LogisticRegression(penalty=penalty,C=c,solver='saga',l1_ratio=.5 if penalty=='elasticnet' else None,max_iter=5000,random_state=3)
  estimators['SVM_'+str(c)]=SVC(C=c,probability=True,random_state=3)
 for depth in [2,4]:
  estimators['RF_'+str(depth)]=RandomForestClassifier(n_estimators=200,max_depth=depth,min_samples_leaf=3,n_jobs=4,random_state=3)
  estimators['GBM_'+str(depth)]=GradientBoostingClassifier(n_estimators=100,max_depth=depth,learning_rate=.05,random_state=3)
 estimators['LDA']=LinearDiscriminantAnalysis(solver='lsqr',shrinkage='auto');estimators['NaiveBayes']=GaussianNB()
 try:
  from xgboost import XGBClassifier
  for depth in [2,4]:estimators['XGBoost_'+str(depth)]=XGBClassifier(n_estimators=150,max_depth=depth,learning_rate=.05,n_jobs=4,random_state=3,eval_metric='logloss')
 except Exception: (O/'XGBoost_error.txt').write_text(traceback.format_exc())
 metrics=[];models={};preds=[];all_oof=[]
 def evaluate(name,model):
  oof=cross_val_predict(model,X,y,cv=cv,method='predict_proba',n_jobs=1)[:,1];fold_auc=[roc_auc_score(y[folds==k],oof[folds==k]) for k in range(10)]
  model.fit(X,y);pr=model.predict_proba(E)[:,1]
  v={'model':name,'CV_mean_fold_AUC':np.mean(fold_auc),'CV_fold_SD':np.std(fold_auc,ddof=1),'CV_pooled_OOF_AUC':roc_auc_score(y,oof),'external_AUC':roc_auc_score(ey,pr),'external_Brier':brier_score_loss(ey,pr)};metrics.append(v);models[name]=model
  for sample,score,target in zip(em.index,pr,ey):preds.append({'model':name,'sample':sample,'y':int(target),'probability':score})
  for sample,score,target,fold in zip(tm.index,oof,y,folds):all_oof.append({'model':name,'sample':sample,'y':int(target),'fold':int(fold),'probability':score})
  pd.DataFrame(metrics).to_csv(O/'model_metrics.tsv',sep='\t',index=False);print('MODEL',v,flush=True)
 for name,est in estimators.items():
  pipe=Pipeline([('scale',StandardScaler()),('RFE',RFE(LogisticRegression(C=1,solver='liblinear',max_iter=2000),n_features_to_select=max(3,len(features)//2),step=1)),('model',est)])
  evaluate(name,pipe)
 # Selection uses training CV only. Internal stacking CV regenerates predictions within each outer training fold.
 top=sorted(metrics,key=lambda z:(-z['CV_mean_fold_AUC'],z['model']))[:3]
 for a,b in itertools.combinations([v['model'] for v in top],2):
  model=StackingClassifier(estimators=[('a',models[a]),('b',models[b])],final_estimator=LogisticRegression(max_iter=2000),cv=StratifiedKFold(5,shuffle=True,random_state=21),n_jobs=1)
  evaluate('Stack_'+a+'__'+b,model)
 best=sorted(metrics,key=lambda z:(-z['CV_mean_fold_AUC'],z['model']))[0];model=models[best['model']]
 pd.DataFrame(preds).to_csv(O/'all_external_predictions.tsv',sep='\t',index=False);pd.DataFrame(all_oof).to_csv(O/'all_training_OOF_predictions.tsv',sep='\t',index=False);joblib.dump(model,O/'chosen_model.joblib')
 (O/'model_selection.json').write_text(json.dumps({'selected':best,'rule':'Highest training CV mean AUC only; external results never used to choose model','limitation':'CV model-selection optimism remains; external cohort is held out. 3 stack pairs selected on global training CV; not a nested assessment of configuration search.'},indent=2))
 import shap
 rng=np.random.default_rng(20261020);background=X[rng.choice(len(X),20,replace=False)]
 explainer=shap.KernelExplainer(lambda v:model.predict_proba(v)[:,1],background)
 values=explainer.shap_values(E,nsamples=100,silent=True)
 np.save(O/'SHAP_all_external.npy',values);pd.DataFrame(values,index=em.index,columns=features).to_csv(O/'SHAP_all_external.tsv',sep='\t')
 pd.DataFrame({'gene':features,'mean_absolute_SHAP':np.abs(values).mean(axis=0)}).sort_values('mean_absolute_SHAP',ascending=False).to_csv(O/'SHAP_importance.tsv',sep='\t',index=False)
 residual=np.max(np.abs(explainer.expected_value+values.sum(axis=1)-model.predict_proba(E)[:,1]));assert residual<1e-5
 v={'status':'completed','train_samples':206,'external_samples':112,'source_TF_features':len(features),'model_configurations':len(metrics),'algorithm_families':9 if any(k.startswith('XGBoost') for k in estimators) else 8,'selected':best,'SHAP_samples':len(E),'SHAP_additivity_max_error':float(residual),'not_executed_original_families':['Stepglm','glmBoost','plsRglm'],'original_exact_113_not_reconstructible':'Unpublished tuning grids and initial feature pool; no fabricated configurations','implementation':'sklearn/XGBoost/SHAP algorithm implementations; original R package versions not claimed'}
 (O/'validation.json').write_text(json.dumps(v,indent=2));print('COMPLETE',json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
