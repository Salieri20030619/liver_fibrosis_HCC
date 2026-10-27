from pathlib import Path
import subprocess,sys,json,traceback
P=Path(__file__).parent;O=P/'F8_F9_real_full_cohort_v1';dep=P/'shap_compatible';dep.mkdir(exist_ok=False)
r=subprocess.run([sys.executable,'-m','pip','install','--target',str(dep),'--no-deps','--disable-pip-version-check','shap'],capture_output=True,text=True,timeout=180);(O/'shap_compat_install.log').write_text(r.stdout+r.stderr);assert r.returncode==0
sys.path.insert(0,str(P/'python_deps'));sys.path.insert(0,str(dep))
try:
 import numpy as np,pandas as pd,joblib,shap
 tr=pd.read_csv(O/'train_features.tsv',sep='\t',index_col=0);ex=pd.read_csv(O/'external_features.tsv',sep='\t',index_col=0)
 metrics=pd.read_csv(O/'model_metrics_with_XGBoost.tsv',sep='\t');best=metrics.sort_values(['CV_mean_fold_AUC','model'],ascending=[False,True]).iloc[0].to_dict()
 model=joblib.load(O/('chosen_model_with_XGBoost.joblib' if str(best['model']).startswith('XGBoost') else 'chosen_model.joblib'))
 X=tr.to_numpy();E=ex.to_numpy();bg=X[np.random.default_rng(20261020).choice(len(X),20,replace=False)]
 explainer=shap.KernelExplainer(lambda a:model.predict_proba(a)[:,1],bg);values=explainer.shap_values(E,nsamples=100,silent=True)
 pd.DataFrame(values,index=ex.index,columns=tr.columns).to_csv(O/'SHAP_all_external.tsv',sep='\t');pd.DataFrame({'gene':tr.columns,'mean_absolute_SHAP':np.abs(values).mean(axis=0)}).sort_values('mean_absolute_SHAP',ascending=False).to_csv(O/'SHAP_importance.tsv',sep='\t',index=False)
 residual=float(np.max(np.abs(explainer.expected_value+values.sum(axis=1)-model.predict_proba(E)[:,1])));assert residual<1e-5
 v={'status':'completed','train_samples':206,'external_samples':112,'source_TF_features':len(tr.columns),'model_configurations':len(metrics),'algorithm_families':9,'selected':best,'SHAP_samples':len(E),'SHAP_version':shap.__version__,'numpy_version':np.__version__,'SHAP_expected_value':float(explainer.expected_value),'SHAP_additivity_max_error':residual,'original_remaining':['Stepglm','glmBoost','plsRglm'],'original_exact113':'Unpublished initial gene pool and grids; reconstructed configurations not exact113','resolved_errors':['Initial import raced dependency install','SHAP0.46 color conversion incompatible with current numpy; installed compatible release in separate project directory, no system mutation']}
 (O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'shap_fix_error.txt').write_text(traceback.format_exc());raise
