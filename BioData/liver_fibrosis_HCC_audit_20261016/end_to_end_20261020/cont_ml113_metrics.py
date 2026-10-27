from pathlib import Path
import json,hashlib,datetime,sys
import numpy as np,pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score,precision_recall_curve,auc,brier_score_loss,confusion_matrix
P=Path(__file__).parent;O=P/'ml113_repaired_20261020_v1';I=O/'inputs';R=O/'results';R.mkdir(exist_ok=False)
terminal=json.loads((O/'workers_terminal.json').read_text());assert terminal['all_returncodes_zero']
methods=pd.read_csv(O/'code/frozen_113_methods.tsv',sep='\t')
tm=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0);em=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0);fold=pd.read_csv(I/'outer_cv_folds.tsv',sep='\t',index_col=0).fold
assert tm.index.equals(fold.index)
full=O/'jobs/full';fv=json.loads((full/'validation.json').read_text());selected=fv['selected_key'];name=fv['selected_method'];threshold=fv['threshold_from_training_OOF']
frozen=json.loads((full/'TRAINING_ONLY_SELECTION_FROZEN.json').read_text())[0];assert frozen['key']==selected
priorbest=pd.read_csv(full/'best_parameters_per_configuration.tsv',sep='\t')
ext=pd.read_csv(full/'all_frozen_external_predictions.tsv.gz',sep='\t',index_col=0);assert ext.index.equals(em.index) and np.array_equal(ext.y,em.y)
outer=pd.DataFrame(np.nan,index=tm.index,columns=methods.method_id);nested=pd.Series(np.nan,index=tm.index);selection=[];fail=[]
for k in range(10):
 q=O/'jobs'/f'outer_{k}';v=json.loads((q/'validation.json').read_text());z=pd.read_csv(q/'all_configuration_outer_predictions.tsv.gz',sep='\t',index_col=0)
 assert z.index.equals(tm.index[fold.eq(k)]) and np.array_equal(z.y,tm.loc[z.index,'y'])
 for key in z.columns.drop('y'):outer.loc[z.index,key.split('|')[0]]=z[key].values
 if v['selected_key'] in z:nested.loc[z.index]=z[v['selected_key']]
 selection.append({'outer_fold':k,'selected_method':v['selected_method'],'selected_key':v['selected_key'],'eligible':v['eligible_configurations'],'heldout_samples':len(z),'all_selected_predictions_finite':np.isfinite(nested.loc[z.index]).all()})
pd.DataFrame(selection).to_csv(R/'nested_outer_selections.tsv',sep='\t',index=False)
outer.to_csv(R/'all113_outer_OOF_probabilities.tsv.gz',sep='\t')
pd.DataFrame({'y':tm.y,'fold':fold,'probability':nested}).to_csv(R/'nested_selection_procedure_OOF.tsv',sep='\t')
def metrics(y,p):
 y=np.asarray(y);p=np.asarray(p)
 if not np.isfinite(p).all():return {'status':'incomplete_predictions','n':len(y),'finite_predictions':int(np.isfinite(p).sum())}
 prec,rec,_=precision_recall_curve(y,p)
 return {'status':'complete','n':len(y),'positive':int(y.sum()),'negative':int((y==0).sum()),'AUROC':float(roc_auc_score(y,p)),'PR_AUC_average_precision':float(average_precision_score(y,p)),'PR_AUC_trapezoidal':float(auc(rec,prec)),'Brier':float(brier_score_loss(y,p)),'constant_same_cohort_prevalence_Brier':float(np.mean((y-y.mean())**2)),'constant_training_prevalence_Brier':float(np.mean((y-tm.y.mean())**2))}
rows=[]
for _,m in methods.iterrows():
 p=outer[m.method_id].to_numpy();r=m.to_dict();r['outer_completed_folds']=int(sum(np.isfinite(p[fold.eq(k)]).all() for k in range(10)))
 r.update({'nested_per_configuration_'+key:value for key,value in metrics(tm.y,p).items()})
 if np.isfinite(p).all():r['nested_per_configuration_mean_fold_AUC']=float(np.mean([roc_auc_score(tm.y[fold.eq(k)],p[fold.eq(k)]) for k in range(10)]))
 z=priorbest[priorbest.method_id.eq(m.method_id)]
 if len(z):
  key=z.iloc[0]['key'];r.update(final_key=key,full_training_CV_mean_AUC=z.iloc[0].mean_fold_AUC,full_training_CV_Brier=z.iloc[0].Brier,selector_parameter=z.iloc[0].selector_parameter,classifier_parameter=z.iloc[0].classifier_parameter)
  if key in ext:r.update({'external_'+key2:value for key2,value in metrics(em.y,ext[key]).items()})
 else:r['status']='No parameter variant completed every full-training CV fold'
 rows.append(r)
allmetrics=pd.DataFrame(rows);allmetrics.to_csv(R/'all113_configuration_results.tsv',sep='\t',index=False)
assert len(allmetrics)==113
failures=[]
for q in (O/'jobs').glob('*/**/variant_status.tsv'):
 if q.relative_to(O/'jobs').parts[0].startswith('smoke'):continue
 z=pd.read_csv(q,sep='\t').fillna('');z=z[z.status.ne('completed')]
 if len(z):
  z['context']=str(q.parent.relative_to(O/'jobs'));failures.append(z)
if failures:
 ff=pd.concat(failures,ignore_index=True);ff.to_csv(R/'all_parameter_failure_records.tsv.gz',sep='\t',index=False)
 ff.groupby(['method_id','configuration','status','error'],dropna=False).size().reset_index(name='failed_partition_parameter_calls').to_csv(R/'failure_summary.tsv',sep='\t',index=False)
else:pd.DataFrame(columns=['method_id','configuration','status','error','failed_partition_parameter_calls']).to_csv(R/'failure_summary.tsv',sep='\t',index=False)
selectedid=selected.split('|')[0];p=ext[selected].to_numpy() if selected in ext else np.full(len(em),np.nan)
assert np.isfinite(p).all(), 'Training-selected final model failed; no external-guided replacement'
trainpred=pd.read_csv(full/'selected_training_OOF.tsv',sep='\t',index_col=0);assert trainpred.index.equals(tm.index)
rng=np.random.default_rng(20261020);pos=np.flatnonzero(em.y.to_numpy()==1);neg=np.flatnonzero(em.y.to_numpy()==0);boot=[];oldboot=[];diff=[]
prior=P/'continuation_newnetwork_models_v2/search';pv=json.loads((prior/'validation.json').read_text());oldname=pv['selected']['model'];old=pd.read_csv(prior/'all243_external_predictions.tsv',sep='\t',index_col=0);assert old.index.equals(em.index);op=old[oldname].to_numpy()
for _ in range(2000):
 ix=np.r_[rng.choice(pos,len(pos),replace=True),rng.choice(neg,len(neg),replace=True)]
 a=roc_auc_score(em.y.to_numpy()[ix],p[ix]);b=roc_auc_score(em.y.to_numpy()[ix],op[ix]);boot.append(a);oldboot.append(b);diff.append(a-b)
fixed=[]
for cohort,y,pp in [('training_selected_variant_CV',tm.y.to_numpy(),trainpred.probability.to_numpy()),('external',em.y.to_numpy(),p)]:
 tn,fp,fn,tp=confusion_matrix(y,pp>=threshold,labels=[0,1]).ravel()
 fixed.append({'cohort':cohort,'threshold':threshold,'TN':int(tn),'FP':int(fp),'FN':int(fn),'TP':int(tp),'sensitivity':tp/(tp+fn),'specificity':tn/(tn+fp),'PPV':tp/(tp+fp) if tp+fp else None,'NPV':tn/(tn+fn) if tn+fn else None,**metrics(y,pp)})
pd.DataFrame(fixed).to_csv(R/'fixed_training_threshold_metrics.tsv',sep='\t',index=False)
pd.DataFrame({'y':em.y,'new113_probability':p,'prior243_probability':op}).to_csv(R/'paired_external_predictions.tsv',sep='\t')
comparison=pd.DataFrame([{'workflow':'repaired113_selection_then_classifier','model':name,'AUC_CI_low':float(np.quantile(boot,.025)),'AUC_CI_high':float(np.quantile(boot,.975)),**metrics(em.y,p)}, {'workflow':'prior243_probability_stacking','model':oldname,'AUC_CI_low':float(np.quantile(oldboot,.025)),'AUC_CI_high':float(np.quantile(oldboot,.975)),**metrics(em.y,op)}]);comparison.to_csv(R/'same_external_cohort_comparison.tsv',sep='\t',index=False)
v={'status':'completed','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selected_method':name,'selected_key':selected,'frozen_selection':frozen,'configurations':113,'complete_nested_configurations':int(allmetrics.outer_completed_folds.eq(10).sum()),'complete_external_configurations':int(allmetrics.get('external_status',pd.Series(dtype=str)).eq('complete').sum()),'nested_selection_procedure':metrics(tm.y,nested),'selected_configuration_nested_OOF':metrics(tm.y,outer[selectedid]),'selected_full_training_tuning_CV':metrics(tm.y,trainpred.probability),'external':metrics(em.y,p),'external_AUC_stratified_bootstrap95CI':np.quantile(boot,[.025,.975]).tolist(),'paired_AUC_difference_new_minus_prior':float(roc_auc_score(em.y,p)-roc_auc_score(em.y,op)),'paired_AUC_difference_bootstrap95CI':np.quantile(diff,[.025,.975]).tolist(),'bootstrap_iterations':2000,'fixed_threshold':threshold,'external_used_for_selection':False,'external_previously_seen_exploratory':True,'limits':'Nested procedure OOF assesses training-only selection among fixed113 and fixed30candidates. Selected-configuration and full-training-tuning CV are selection-conditional; external was previously seen and motivated adaptation. Cross-sectional F3/F4, not incidentHCC. No formalfigures.'}
(R/'validation.json').write_text(json.dumps(v,indent=2,allow_nan=False));print(json.dumps(v),flush=True)
