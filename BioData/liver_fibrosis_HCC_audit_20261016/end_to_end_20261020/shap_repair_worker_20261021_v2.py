from pathlib import Path
import sys,os,json,pickle,hashlib,time,traceback
P=Path(__file__).parent;O=P/'SHAP_repair_243_20261021_v1';sys.path[:0]=[str(O/'code'),str(P/'shap_compatible'),str(P/'python_deps')]
import numpy as np,pandas as pd,joblib,shap
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from frozen_explainer import FrozenStack,kernel_explainer,explain_row
from test_frozen_explainer import run_tests
N=P/'continuation_newnetwork_models_v2';I=N/'features'
def js(d,name):(O/name).write_text(json.dumps(d,indent=2,default=str))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
try:
 X=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0);E=pd.read_csv(I/'external_features.tsv',sep='\t',index_col=0)
 models={n:joblib.load(I/(n+'.joblib')) for n in ['RF_4','SVM_1.0']}
 with open(N/'search/meta_models.pkl','rb') as f:meta=pickle.load(f)['Stack_RF_4__SVM_1.0']
 stack=FrozenStack(models,meta,X);bg=X.sample(20,random_state=21)
 bg.to_csv(O/'frozen_training_background.tsv',sep='\t')
 files=[I/'train_features.tsv',I/'external_features.tsv',I/'train_design.tsv',I/'external_design.tsv',I/'outer_cv_folds.tsv',I/'RF_4.joblib',I/'SVM_1.0.joblib',N/'search/meta_models.pkl',N/'search/all243_external_predictions.tsv',N/'SHAP/external_gene_SHAP.tsv']
 manifest=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in files];pd.DataFrame(manifest).to_csv(O/'frozen_source_manifest.tsv',sep='\t',index=False)
 rule={'scale':'positive-class probability; identity link','algorithm':'KernelSHAP on exact15-dimensional dummy-feature reduced interventional game; unregularized weighted least squares l1_reg=0; lift by dummy-player identity to original30 columns','background':'same20 training rows as historical output, pandas sample random_state21, uniform weights; no external labels','budgets':[512,1024,2048,4096],'seeds':[20261021,20261022],'precision_acceptance':'per cohort active-feature importance rank Spearman>=.98, top5 overlap>=4, relative importance L1 difference<=.05, sample-wise attribution MAE<=.002 probability; both budget doubling and seed comparison','max_budget':4096,'active':stack.active,'full_columns':stack.columns,'shap_version':shap.__version__,'training_attribution':'final frozen model on training rows; in-sample explanation, not OOF explanation','no_new_model_or_feature_selection':True}
 js(rule,'frozen_explanation_rules.json')
 saved=pd.read_csv(N/'search/all243_external_predictions.tsv',sep='\t',index_col=0).loc[E.index,'Stack_RF_4__SVM_1.0'].to_numpy()
 with open(O/'regression_tests.txt','w') as log:tests=run_tests(stack,X,E,saved,log)
 assert tests.wasSuccessful()
 js({'tests_run':tests.testsRun,'failures':len(tests.failures),'errors':len(tests.errors),'passed':tests.wasSuccessful()},'regression_tests.json')
 exp=kernel_explainer(stack,bg);baseline=float(exp.expected_value)
 predictions={c:stack.predict_full(x) for c,x in [('training',X),('external',E)]}
 for c,x in [('training',X),('external',E)]:pd.DataFrame({'sample':x.index,'probability':predictions[c]}).to_csv(O/(c+'_frozen_predictions.tsv'),sep='\t',index=False)
 tm=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0);em=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0)
 js({'expected_value':baseline,'external_prediction_max_error':float(np.max(abs(predictions['external']-saved))),'external_AUC_unchanged':float(roc_auc_score(em.loc[E.index,'y'],predictions['external'])),'train_model_fitted_AUC_not_CV':float(roc_auc_score(tm.loc[X.index,'y'],predictions['training']))},'prediction_validation.json')
 records=[];previous={};final={};start=time.time()
 def compare(a,b,cohort,kind,budget):
  av=abs(a.loc[:,stack.active]).mean();bv=abs(b.loc[:,stack.active]).mean();rho=float(spearmanr(av,bv).statistic);overlap=len(set(av.nlargest(5).index)&set(bv.nlargest(5).index));rel=float(abs(av-bv).sum()/max(bv.sum(),1e-12));mae=float(abs(a.loc[:,stack.active]-b.loc[:,stack.active]).to_numpy().mean());maxerr=float(abs(a-b).to_numpy().max())
  d=dict(cohort=cohort,comparison=kind,budget=budget,rank_spearman=rho,top5_overlap=overlap,relative_importance_L1=rel,sample_gene_MAE=mae,max_sample_gene_difference=maxerr,passed=rho>=.98 and overlap>=4 and rel<=.05 and mae<=.002);records.append(d);return d['passed']
 for budget in rule['budgets']:
  current={}
  for seed in rule['seeds']:
   for cohort,x in [('training',X),('external',E)]:
    path=O/(f'{cohort}_budget{budget}_seed{seed}.tsv');values=[]
    for j in range(len(x)):
     row=explain_row(exp,stack,x.iloc[j:j+1],budget,seed+j*997);values.append(row)
     if (j+1)%25==0:
      js({'stage':'explaining','budget':budget,'seed':seed,'cohort':cohort,'completed_samples':j+1,'total_samples':len(x),'elapsed_seconds':time.time()-start},'progress.json');print(budget,seed,cohort,j+1,'elapsed',round(time.time()-start),flush=True)
    v=pd.concat(values);assert v.index.equals(x.index) and list(v.columns)==list(X.columns)
    assert abs(v.loc[:,~stack.mask]).to_numpy().max()==0
    assert np.max(abs(baseline+v.sum(1).to_numpy()-predictions[cohort]))<1e-8
    v.to_csv(path,sep='\t');current[(seed,cohort)]=v
  ok=[]
  for cohort in ['training','external']:
   ok.append(compare(current[(20261021,cohort)],current[(20261022,cohort)],cohort,'independent_seeds',budget))
   if previous:ok.append(compare(previous[(20261021,cohort)],current[(20261021,cohort)],cohort,'budget_doubling',budget))
  pd.DataFrame(records).to_csv(O/'convergence_checks.tsv',sep='\t',index=False)
  if previous and all(ok):final=current;break
  previous=current;final=current
 convergence=bool(previous and all(ok) and budget>512)
 # Average independent MC solutions at the accepted budget; preserve each replicate.
 summaries=[];errs={}
 for cohort,x in [('training',X),('external',E)]:
  v=(final[(20261021,cohort)]+final[(20261022,cohort)])/2;v.to_csv(O/(cohort+'_repaired30_SHAP.tsv'),sep='\t')
  importance=abs(v).mean();out=pd.DataFrame({'gene':X.columns,'used_by_model':stack.mask,'mean_abs_SHAP':importance.values,'mean_signed_SHAP':v.mean().values,'rank':importance.rank(ascending=False,method='min').astype(int).values,'cohort':cohort})
  out.to_csv(O/(cohort+'_importance.tsv'),sep='\t',index=False);summaries.append(out)
  er=baseline+v.sum(1).to_numpy()-predictions[cohort];errs[cohort]=float(np.max(abs(er)));pd.DataFrame({'sample':x.index,'expected_value':baseline,'SHAP_sum':v.sum(1).to_numpy(),'frozen_probability':predictions[cohort],'reconstruction_error':er}).to_csv(O/(cohort+'_local_reconstruction.tsv'),sep='\t',index=False)
 old=pd.read_csv(N/'SHAP/gene_importance.tsv',sep='\t',index_col=0).iloc[:,0];new=summaries[1].set_index('gene');comparison=new.copy();comparison['old_mean_abs_SHAP']=old;comparison['old_rank']=old.rank(ascending=False,method='min').astype(int);comparison['rank_change_new_minus_old']=comparison['rank']-comparison.old_rank;comparison.reset_index().to_csv(O/'old_new_external_rank_comparison.tsv',sep='\t',index=False)
 assert all(sha(Path(r['path']))==r['sha256'] for r in manifest)
 js({'status':'completed' if convergence else 'precision_not_converged_at_budget_cap','convergence_passed':convergence,'final_budget_per_seed':budget,'seeds':rule['seeds'],'combined_by':'mean of two independent seed attributions, NOT average of base-model SHAP','tests_passed':True,'prediction_error':float(np.max(abs(predictions['external']-saved))),'reconstruction_max_error':errs,'dummy_features_zero':True,'dimensions':[206,112,30],'source_hashes_unchanged':True,'elapsed_seconds':time.time()-start,'external_top5':new.mean_abs_SHAP.nlargest(5).index.tolist(),'training_top5':summaries[0].nlargest(5,'mean_abs_SHAP').gene.tolist()},'validation.json')
 js({'historical_output':str(N/'SHAP'),'historical_files_preserved':True,'status':'superseded_for_attribution_interpretation' if convergence else 'historical_invalid_dummy_property_new_precision_pending','replacement':str(O),'reason':'Historical30-dimensional100-sample KernelSHAP with all30 regression terms assigned nonzero effects to15 functionally absent inputs; finite coalition approximation plus surplus regressors. New game reduces exact dummy players before estimation, disables feature-selection regularization, validates two seeds/budgets.'},'SUPERSESSION.json')
 print((O/'validation.json').read_text(),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
