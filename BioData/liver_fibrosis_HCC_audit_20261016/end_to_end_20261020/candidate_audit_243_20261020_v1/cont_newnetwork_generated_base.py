from pathlib import Path
import sys,json,itertools,traceback
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));O=P/'continuation_newnetwork_models_v2'/'search';O.mkdir(exist_ok=True);I=P/'continuation_newnetwork_models_v2'/'features'
import numpy as np,pandas as pd
from sklearn.model_selection import StratifiedKFold,cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from xgboost import XGBClassifier
X=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0);y=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0).loc[X.index,'y'].to_numpy();fold=pd.read_csv(I/'outer_cv_folds.tsv',sep='\t').set_index('sample').loc[X.index,'fold'].to_numpy();est={}
for c in [.1,1.]:
 for name,penalty in [('Lasso','l1'),('Ridge','l2'),('ElasticNet','elasticnet')]:est[name+'_'+str(c)]=LogisticRegression(penalty=penalty,C=c,solver='saga',l1_ratio=.5 if penalty=='elasticnet' else None,max_iter=5000,random_state=3)
 est['SVM_'+str(c)]=SVC(C=c,probability=True,random_state=3)
for depth in [2,4]:
 est['RF_'+str(depth)]=RandomForestClassifier(n_estimators=200,max_depth=depth,min_samples_leaf=3,n_jobs=4,random_state=3)
 est['GBM_'+str(depth)]=GradientBoostingClassifier(n_estimators=100,max_depth=depth,learning_rate=.05,random_state=3)
 est['XGBoost_'+str(depth)]=XGBClassifier(n_estimators=150,max_depth=depth,learning_rate=.05,n_jobs=4,random_state=3,eval_metric='logloss')
est['LDA']=LinearDiscriminantAnalysis(solver='lsqr',shrinkage='auto');est['NaiveBayes']=GaussianNB()

import joblib
from sklearn.metrics import roc_auc_score,brier_score_loss
E=pd.read_csv(I/'external_features.tsv',sep='\t',index_col=0);ey=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0).loc[E.index,'y'].to_numpy();metrics=[];trrows=[];erows=[]
cv=[(np.flatnonzero(fold!=k),np.flatnonzero(fold==k)) for k in range(10)]
for name,model in est.items():
 pipe=Pipeline([('scale',StandardScaler()),('RFE',RFE(LogisticRegression(C=1,solver='liblinear',max_iter=2000),n_features_to_select=max(3,X.shape[1]//2),step=1)),('model',model)])
 p=cross_val_predict(pipe,X.to_numpy(),y,cv=cv,method='predict_proba',n_jobs=1)[:,1];pipe.fit(X.to_numpy(),y);q=pipe.predict_proba(E.to_numpy())[:,1];joblib.dump(pipe,I/(name+'.joblib'));auc=[roc_auc_score(y[fold==k],p[fold==k]) for k in range(10)];metrics.append({'model':name,'CV_mean_fold_AUC':np.mean(auc),'external_AUC':roc_auc_score(ey,q),'external_Brier':brier_score_loss(ey,q)})
 trrows.extend(dict(model=name,sample=s,y=int(yy),fold=int(k),probability=float(v)) for s,yy,k,v in zip(X.index,y,fold,p));erows.extend(dict(model=name,sample=s,y=int(yy),probability=float(v)) for s,yy,v in zip(E.index,ey,q));print('base',name,flush=True)
for records,stem in [(trrows,'training_OOF'),(erows,'external')]:
 d=pd.DataFrame(records);d[~d.model.str.startswith('XGBoost')].to_csv(I/('all_'+stem+'_predictions.tsv'),sep='\t',index=False);d[d.model.str.startswith('XGBoost')].to_csv(I/('XGBoost_'+stem+'_predictions.tsv'),sep='\t',index=False)
pd.DataFrame(metrics).to_csv(I/'base_model_metrics.tsv',sep='\t',index=False);(I/'PYTHON_BASE_COMPLETE.txt').write_text('completed')
