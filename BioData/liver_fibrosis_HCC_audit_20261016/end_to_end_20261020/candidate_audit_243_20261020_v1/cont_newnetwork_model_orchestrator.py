from pathlib import Path
import os,time,json,subprocess,hashlib,shutil
import numpy as np,pandas as pd
P=Path(__file__).parent;N=P/'continuation_newnetwork_models_v2';N.mkdir(exist_ok=False);I=N/'features';I.mkdir();C=P/'continuation_SCENIC_candidates_v1';old=P/'F8_F9_real_full_cohort_v1'
for i in range(1080):
 if (C/'validation.json').exists():break
 if (P/'continuation_SCENIC_expanded_v1/error.txt').exists():raise RuntimeError('Newnetwork upstream failed, no downstreammodels attempted')
 time.sleep(10)
else:raise RuntimeError('Newnetwork candidates pending')
tf=set(pd.read_csv(C/'F5G_candidates.tsv',sep='\t').TF);F=P/'F6_S4_S6_F8_gene_validation';tr=pd.read_csv(F/'GSE135251_log2CPM_full.tsv.gz',sep='\t',index_col=0);ex=pd.read_csv(F/'GSE162694_log2CPM_full.tsv.gz',sep='\t',index_col=0);tm=pd.read_csv(old/'train_design.tsv',sep='\t',index_col=0);em=pd.read_csv(old/'external_design.tsv',sep='\t',index_col=0);tr=tr[tm.index];ex=ex[em.index];genes=sorted(tf&set(tr.index)&set(ex.index));genes=[g for g in genes if tr.loc[g].std()>0];prov={'source_TF_candidate_sha256':hashlib.sha256((C/'F5G_candidates.tsv').read_bytes()).hexdigest(),'new_network_TFs':sorted(tf),'matched_nonconstant_features':genes,'unmatched_or_constant':sorted(tf-set(genes)),'endpoint':'crosssectionalF3/F4versusF0/F1/F2,notincidentHCC','external_status':'Previouslyseenexternalcohort; exploratoryexpandedvalidation,notpristineconfirmatory'};(I/'feature_provenance.json').write_text(json.dumps(prov,indent=2))
if len(genes)<3:
 (N/'validation.json').write_text(json.dumps({'status':'not_estimable_fewer_than3_available_TFs','available_features':genes},indent=2));raise SystemExit(0)
X=tr.loc[genes].T;E=ex.loc[genes].T;X.to_csv(I/'train_features.tsv',sep='\t');E.to_csv(I/'external_features.tsv',sep='\t');tm.to_csv(I/'train_design.tsv',sep='\t');em.to_csv(I/'external_design.tsv',sep='\t');shutil.copyfile(old/'outer_cv_folds.tsv',I/'outer_cv_folds.tsv');prior=pd.read_csv(old/'train_features.tsv',sep='\t',index_col=0)
if genes==prior.columns.tolist() and np.allclose(X.to_numpy(),prior.loc[X.index].to_numpy(),atol=1e-10):
 (N/'validation.json').write_text(json.dumps({'status':'reused_identical_feature_matrix','source':str(P/'continuation_stack_search_v1'),'reason':'Newnetworkavailablecandidateinterface identical; reuse completed243models andSHAP without duplicatejobs'},indent=2));raise SystemExit(0)
bundle=json.loads((P/'newnetwork_model_script_bundle.json').read_text());paths={}
for name,source in bundle.items():
 q=P/('cont_newnetwork_generated_'+name+'.py');assert not q.exists();q.write_text(source);paths[name]=q
env=os.environ.copy();env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4',NUMBA_NUM_THREADS='4');jobs={};logs=[]
for name in ['base','python_inner','R_base','R_inner','combine']:
 f=open(N/(name+'.log'),'w');logs.append(f);job=subprocess.Popen([os.sys.executable,str(paths[name])],stdout=f,stderr=subprocess.STDOUT,env=env,cwd=P);jobs[name]=job
(N/'jobs.json').write_text(json.dumps({name:job.pid for name,job in jobs.items()},indent=2));print('Newfeatures',len(genes),'launched',list(jobs),flush=True)
returns={name:job.wait() for name,job in jobs.items()}
if any(returns.values()):
 (N/'validation.json').write_text(json.dumps({'status':'failed_component','exit_codes':returns},indent=2));raise RuntimeError(str(returns))
with open(N/'metrics.log','w') as f:z=subprocess.run([os.sys.executable,str(paths['metrics'])],env=env,stdout=f,stderr=subprocess.STDOUT)
v={'status':'models_and_metrics_complete_SHAP_pending','feature_count':len(genes),'configurations':243,'exit_codes':returns,'metric_returncode':z.returncode,'source_candidate_file':str(C/'F5G_candidates.tsv'),'changed_feature_pool':True};(N/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
