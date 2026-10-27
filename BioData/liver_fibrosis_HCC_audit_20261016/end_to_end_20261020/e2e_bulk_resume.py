from pathlib import Path
import sys,json,traceback
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));O=P/'F8_F9_real_full_cohort_v1'
import numpy as np,pandas as pd,joblib
from sklearn.model_selection import StratifiedKFold,cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,brier_score_loss
from xgboost import XGBClassifier
try:
 tr=pd.read_csv(O/'train_features.tsv',sep='\t',index_col=0);ex=pd.read_csv(O/'external_features.tsv',sep='\t',index_col=0);tm=pd.read_csv(O/'train_design.tsv',sep='\t',index_col=0);em=pd.read_csv(O/'external_design.tsv',sep='\t',index_col=0);X=tr.to_numpy();E=ex.to_numpy();y=tm.y.to_numpy();ey=em.y.to_numpy();features=tr.columns.tolist()
 metrics=pd.read_csv(O/'model_metrics.tsv',sep='\t').to_dict('records');models={};cv=StratifiedKFold(10,shuffle=True,random_state=20261020);folds=pd.read_csv(O/'outer_cv_folds.tsv',sep='\t').fold.to_numpy();preds=[];oofs=[]
 for depth in [2,4]:
  name='XGBoost_'+str(depth);est=XGBClassifier(n_estimators=150,max_depth=depth,learning_rate=.05,n_jobs=4,random_state=3,eval_metric='logloss')
  pipe=Pipeline([('scale',StandardScaler()),('RFE',RFE(LogisticRegression(C=1,solver='liblinear',max_iter=2000),n_features_to_select=max(3,len(features)//2))),('model',est)])
  oof=cross_val_predict(pipe,X,y,cv=cv,method='predict_proba')[:,1];auc=[roc_auc_score(y[folds==k],oof[folds==k]) for k in range(10)];pipe.fit(X,y);pr=pipe.predict_proba(E)[:,1];models[name]=pipe
  v={'model':name,'CV_mean_fold_AUC':float(np.mean(auc)),'CV_fold_SD':float(np.std(auc,ddof=1)),'CV_pooled_OOF_AUC':roc_auc_score(y,oof),'external_AUC':roc_auc_score(ey,pr),'external_Brier':brier_score_loss(ey,pr)};metrics.append(v);print(v,flush=True)
  preds.extend({'model':name,'sample':s,'y':int(t),'probability':float(v)} for s,t,v in zip(em.index,ey,pr));oofs.extend({'model':name,'sample':s,'y':int(t),'fold':int(f),'probability':float(v)} for s,t,f,v in zip(tm.index,y,folds,oof))
 pd.DataFrame(metrics).to_csv(O/'model_metrics_with_XGBoost.tsv',sep='\t',index=False)
 pd.DataFrame(preds).to_csv(O/'XGBoost_external_predictions.tsv',sep='\t',index=False);pd.DataFrame(oofs).to_csv(O/'XGBoost_training_OOF_predictions.tsv',sep='\t',index=False)
 best=sorted(metrics,key=lambda z:(-z['CV_mean_fold_AUC'],z['model']))[0];model=models.get(best['model'])
 if model is None:model=joblib.load(O/'chosen_model.joblib')
 else:joblib.dump(model,O/'chosen_model_with_XGBoost.joblib')
 import shap
 bg=X[np.random.default_rng(20261020).choice(len(X),20,replace=False)];explainer=shap.KernelExplainer(lambda a:model.predict_proba(a)[:,1],bg)
 values=explainer.shap_values(E,nsamples=100,silent=True)
 pd.DataFrame(values,index=ex.index,columns=features).to_csv(O/'SHAP_all_external.tsv',sep='\t');pd.DataFrame({'gene':features,'mean_absolute_SHAP':np.abs(values).mean(axis=0)}).sort_values('mean_absolute_SHAP',ascending=False).to_csv(O/'SHAP_importance.tsv',sep='\t',index=False)
 residual=float(np.max(np.abs(explainer.expected_value+values.sum(axis=1)-model.predict_proba(E)[:,1])));assert residual<1e-5
 v={'status':'completed','train_samples':206,'external_samples':112,'source_TF_features':len(features),'model_configurations':len(metrics),'algorithm_families':9,'selected':best,'SHAP_samples':len(E),'SHAP_expected_value':float(explainer.expected_value),'SHAP_additivity_max_error':residual,'original_remaining':['Stepglm','glmBoost','plsRglm'],'original_exact113':'Unpublished initial gene pool and grids; reconstructed19 configurations not exact113','resolved_error':'Initial SHAP import raced dependency installation; resumed from saved trained model, did not discard results.'}
 (O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'resume_error.txt').write_text(traceback.format_exc());raise
