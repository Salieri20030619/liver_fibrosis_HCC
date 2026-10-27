from pathlib import Path
import sys,json,time,itertools,pickle
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'));O=P/'continuation_newnetwork_models_v2'/'search';I=P/'continuation_newnetwork_models_v2'/'features';R=P/'continuation_newnetwork_models_v2'/'R_base'
import numpy as np,pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,brier_score_loss,average_precision_score
for attempt in range(180):
 if (O/'frozen_search_rules.json').exists():break
 time.sleep(2)
rules=json.loads((O/'frozen_search_rules.json').read_text());models=rules['base_models'];pairs=rules['pairs']
for retry in range(360):
 if (O/'R_cache_COMPLETE.txt').exists() and (O/'python_cache_complete.json').exists() and (I/'PYTHON_BASE_COMPLETE.txt').exists() and (R/'validation.json').exists():break
 time.sleep(10)
else:raise RuntimeError('Caches incomplete; combination search not executed')
tr=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0);ex=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0);fold=pd.read_csv(I/'outer_cv_folds.tsv',sep='\t').set_index('sample').loc[tr.index,'fold'];y=tr.y.to_numpy();ey=ex.y.to_numpy()
def load(files,index):
 a=pd.concat([pd.read_csv(f,sep='\t') for f in files]);a=a[a.model.isin(models)];assert not a.duplicated(['sample','model']).any();a=a.pivot(index='sample',columns='model',values='probability').loc[index,models];assert np.isfinite(a.to_numpy()).all();return a
base=load([I/'all_training_OOF_predictions.tsv',I/'XGBoost_training_OOF_predictions.tsv',R/'training_OOF_predictions.tsv'],tr.index)
ext=load([I/'all_external_predictions.tsv',I/'XGBoost_external_predictions.tsv',R/'external_predictions.tsv'],ex.index)
cache={}
for c in ['full']+list(range(10)):
 a=pd.read_csv(O/f'python_inner_oof_{c}.tsv',sep='\t',index_col=0);b=pd.read_csv(O/f'R_inner_oof_{c}.tsv',sep='\t',index_col=0);cache[c]=a.join(b).loc[:,models];assert np.isfinite(cache[c]).all().all()
metrics=[];cvpred={};epred={};metas={}
def record(name,p,q,kind):
 aucs=[roc_auc_score(y[fold==k],p[fold==k]) for k in range(10)];metrics.append(dict(model=name,kind=kind,CV_mean_fold_AUC=np.mean(aucs),CV_fold_SD=np.std(aucs,ddof=1),CV_pooled_OOF_AUC=roc_auc_score(y,p),external_AUC=roc_auc_score(ey,q),external_Brier=brier_score_loss(ey,q),external_AP=average_precision_score(ey,q)));cvpred[name]=p;epred[name]=q
for name in models:record(name,base[name].to_numpy(),ext[name].to_numpy(),'single')
for a,b in pairs:
 name='Stack_'+a+'__'+b;p=np.zeros(len(y))
 for k in range(10):
  ids=tr.index[fold!=k];test=fold==k;fit=LogisticRegression(max_iter=2000).fit(cache[k].loc[ids,[a,b]],tr.loc[ids,'y']);p[test]=fit.predict_proba(base.loc[test,[a,b]])[:,1]
 fit=LogisticRegression(max_iter=2000).fit(cache['full'].loc[tr.index,[a,b]],y);q=fit.predict_proba(ext[[a,b]])[:,1];metas[name]=fit;record(name,p,q,'stack')
m=pd.DataFrame(metrics).sort_values(['CV_mean_fold_AUC','model'],ascending=[False,True]);assert len(m)==243
m.to_csv(O/'all243_model_metrics.tsv',sep='\t',index=False);pd.DataFrame(cvpred,index=tr.index).to_csv(O/'all243_training_OOF.tsv',sep='\t');pd.DataFrame(epred,index=ex.index).to_csv(O/'all243_external_predictions.tsv',sep='\t')
best=m.iloc[0].model;bestq=epred[best];rng=np.random.default_rng(20261020);aucs=[]
for i in range(2000):
 ids=rng.integers(0,len(ey),len(ey))
 if len(np.unique(ey[ids]))==2:aucs.append(roc_auc_score(ey[ids],bestq[ids]))
with open(O/'meta_models.pkl','wb') as f:pickle.dump(metas,f)
v={'status':'completed','configurations':243,'single':22,'cross_family_pairs':221,'selected_by':'training_mean_outer_fold_AUC_only','selected':m.iloc[0].to_dict(),'external_AUC_bootstrap_CI':np.quantile(aucs,[.025,.975]).tolist(),'external_constant_prevalence_Brier':float(np.var(ey)),'important_limits':rules['original113_difference']+' '+rules['validation'],'upstream':'NewnetworkTF candidates from actual expanded original-direction RSS intersection; actual previously fitted outer/base predictions reused; all combination meta learners trained on newly computed inner5 OOF.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
