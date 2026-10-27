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
models=list(est)+['Stepglm_backward','Stepglm_both','glmBoost_100','glmBoost_200','plsRglm_2','plsRglm_3'];pairs=[(a,b) for a,b in itertools.combinations(sorted(models),2) if a.split('_')[0]!=b.split('_')[0]];assert len(models)==22 and len(pairs)==221
rules={'status':'frozen_before_new_search','feature_count':X.shape[1],'base_models':models,'pairs':pairs,'total_configurations':243,'pair_rule':'All unordered pairs of existing22 frozen single variants from different algorithm families; excludes10 within-family pairs. Not a fabricated91 subset.','original113_difference':'Original91 exact pair identities/grid not supplied in available methods. Transparent complete reconstruction=22+221, not exact113.','validation':'Outer10 folds unchanged; meta learner trained solely on inner5 OOF within each outer training subset. Fulltraining inner5 for final fit; externalneverselects. All22 candidatepool comesfromnewnetwork persistedTFintersection; RFE/scaling nested. Search optimism still exists for selecting among243 viaCV; externalalreadyexamined in previousphase so not pristine confirmatory evaluation.','identity':'Bulk206 source sample identifiers; external112 expressionIDs, person-level crosssource independence not proven.'};(O/'frozen_search_rules.json').write_text(json.dumps(rules,indent=2))
allrows=[]
for outer in ['full']+list(range(10)):
 use=np.ones(len(y),bool) if outer=='full' else fold!=outer;ids=np.flatnonzero(use);cv=StratifiedKFold(5,shuffle=True,random_state=21);out=pd.DataFrame(index=X.index[use]);out.index.name='sample'
 # Publish splits so R variants use exactly the same inner folds.
 inner=np.zeros(len(ids),int)
 for k,(_,test) in enumerate(cv.split(X.to_numpy()[use],y[use])):inner[test]=k
 pd.DataFrame({'sample':out.index,'inner_fold':inner}).to_csv(O/('inner_folds_'+str(outer)+'.tsv'),sep='\t',index=False)
 for name,model in est.items():
  pipe=Pipeline([('scale',StandardScaler()),('RFE',RFE(LogisticRegression(C=1,solver='liblinear',max_iter=2000),n_features_to_select=max(3,X.shape[1]//2),step=1)),('model',model)])
  out[name]=cross_val_predict(pipe,X.to_numpy()[use],y[use],cv=cv,method='predict_proba',n_jobs=1)[:,1];assert np.isfinite(out[name]).all();out.to_csv(O/('python_inner_oof_'+str(outer)+'.tsv'),sep='\t');print('inner cache',outer,name,flush=True)
 (O/('python_'+str(outer)+'_COMPLETE.txt')).write_text('complete')
(O/'python_cache_complete.json').write_text(json.dumps({'status':'completed','outer_contexts':11,'base_models':16},indent=2))
