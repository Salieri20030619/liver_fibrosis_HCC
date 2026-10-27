from pathlib import Path
import json,sys,hashlib,gzip,tarfile,datetime,pickle
import numpy as np,pandas as pd
from scipy.stats import spearmanr,mannwhitneyu
P=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/end_to_end_20261020');D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'candidate_audit_243_20261020_v1';N=P/'continuation_newnetwork_models_v2';I=N/'features'
assert json.loads((O/'validation.json').read_text())['status']=='completed'
sys.path.insert(0,str(P/'python_deps'));import joblib
X=pd.read_csv(I/'train_features.tsv',sep='\t',index_col=0);E=pd.read_csv(I/'external_features.tsv',sep='\t',index_col=0);tm=pd.read_csv(I/'train_design.tsv',sep='\t',index_col=0);em=pd.read_csv(I/'external_design.tsv',sep='\t',index_col=0)
def save(d,n):d.to_csv(O/n,sep='\t',index=False)
def js(d,n):(O/n).write_text(json.dumps(d,indent=2,default=str))
# Preserve public source protocol/identity fields only, no contact or identity guessing.
records=json.loads((P/'GSE135251_parsed_samples.json').read_text())
fields=['GSM','title','source_name_ch1','characteristics_ch1','description','data_processing','platform_id','relation','molecule_ch1','extract_protocol_ch1']
js([{k:r.get(k,[]) for k in fields} for r in records if r['GSM'] in X.index],'train_public_metadata.json')
m=pd.read_csv(P/'continuation_network_communication_inputs_v1/SCENIC/cell_metadata.tsv',sep='\t',index_col=0)
sets={'network':set(m.GSM.dropna()),'train':set(tm.index),'external':set(em.GSM)}
js({'GSM_counts':{k:len(v) for k,v in sets.items()},'GSM_overlaps':{a+'__'+b:sorted(sets[a]&sets[b]) for a,b in [('network','train'),('network','external'),('train','external')]},'train_donor_label_unique':int(tm.donor.nunique()),'external_expression_ID_unique':int(E.index.nunique()),'external_GSM_unique':int(em.GSM.nunique()),'train_duplicate_feature_vectors':int(X.duplicated().sum()),'external_duplicate_feature_vectors':int(E.duplicated().sum()),'limitations':'No shared GSM is not proof of no shared people. Patient IDs are study-scoped and cannot be equated by removing P/patient prefixes. No validated genotype/author cross-study identity key.'},'sample_overlap_checks.json')
# Saved prediction domain diagnostics.
pp=[]
for name in ['RF_4','SVM_1.0']:
 mod=joblib.load(I/(name+'.joblib'));sc=mod.named_steps['scale'];mask=mod.named_steps['RFE'].support_
 for label,x in [('train_fitted',X),('external',E)]:
  q=mod.predict_proba(x.to_numpy())[:,1];z=sc.transform(x.to_numpy())[:,mask]
  pp.append(dict(model=name,cohort=label,prob_min=float(q.min()),prob_q05=float(np.quantile(q,.05)),prob_median=float(np.median(q)),prob_q95=float(np.quantile(q,.95)),prob_max=float(q.max()),mean_abs_selected_z=float(np.mean(abs(z))),fraction_selected_z_abs_gt3=float(np.mean(abs(z)>3))))
save(pd.DataFrame(pp),'saved_model_prediction_domain.tsv')
sc=joblib.load(I/'RF_4.joblib').named_steps['scale'];ez=pd.DataFrame(sc.transform(E.to_numpy()),columns=E.columns,index=E.index)
save(pd.DataFrame({'gene':X.columns,'external_mean_training_z':ez.mean().values,'external_SD_training_z':ez.std().values,'external_fraction_abs_z_gt3':(abs(ez)>3).mean().values}),'feature_standardization_audit.tsv')
sex=[]
for s,g in em.groupby('sex'):
 for gene in X.columns:
  v=E.loc[g.index,gene];sex.append(dict(gene=gene,sex=s,n=len(g),advanced=int(g.y.sum()),mean=float(v.mean()),median=float(v.median())))
save(pd.DataFrame(sex),'external_sex_feature_diagnostics.tsv')
stage=em.groupby(['sex','F']).size().rename('n').reset_index();save(stage,'external_sex_stage_counts.tsv')
anno=P/'continuation_annotation_refined_v1/subcluster_complete_marker_evidence.tsv.gz'
if anno.exists():
 markers=pd.read_csv(anno,sep='\t');save(markers[markers.parent_lineage.eq('Myeloid')],'myeloid_state_marker_evidence.tsv')
diag=P/'continuation_annotation_evidence_v1/cell_lineage_codetection_diagnostics.tsv.gz'
if diag.exists():
 z=pd.read_csv(diag,sep='\t',index_col=0).loc[m.index];cols=[c for c in z if ('fraction' in c or 'flag' in c or 'markers_detected' in c)];z['state']=m.refined_cluster;z['study']=m.study
 save(z.groupby(['state','study'])[cols].mean().reset_index(),'sampled_myeloid_lineage_codetection.tsv')
# Record contribution concentration in the advanced discovery cells.
q=m[m.regulon_role.eq('discovery_advanced_F3_F4')];out=[]
for state,g in q.groupby('refined_cluster'):
 sizes=g.source_record.value_counts();w=sizes/sizes.sum()
 out.append(dict(state=state,cells=len(g),source_records=len(sizes),min_record_cells=int(sizes.min()),max_record_cells=int(sizes.max()),largest_record_fraction=float(w.max()),effective_equal_cell_contribution_records=float(1/(w*w).sum()),study_counts=json.dumps(g.study.value_counts().to_dict()),F_counts=json.dumps(g.F.value_counts().to_dict())))
save(pd.DataFrame(out),'discovery_state_contribution_concentration.tsv')
# Exact reference files and runtime defaults from the installed CLI help.
import subprocess
V=D/'F5_setup_v1/pyscenic_venv/bin/pyscenic'
for stage in ['grn','ctx','aucell']:
 r=subprocess.run([str(V),stage,'--help'],capture_output=True,text=True)
 (O/('pyscenic_'+stage+'_help.txt')).write_text(r.stdout+r.stderr)
ref=[D/'F5_regulatory_resource_audit_v1/allTFs_hg38.txt']+list((P.parent/'reference_cisTarget_hg38_mc9nr_v3').glob('*.feather'))+list((P.parent/'reference_cisTarget_aux_v1').glob('*.tbl'))
save(pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in ref]),'reference_file_manifest.tsv')
for p in [P/'cand_audit_compute_v1.py',P/'cont_SCENIC_expanded.py',P/'cont_SCENIC_postprocess.py',P/'cont_network_inputs_review_gated.py',P/'cont_newnetwork_model_orchestrator.py',P/'cont_newnetwork_generated_base.py',P/'cont_newnetwork_generated_combine.py',P/'cont_newnetwork_selected_SHAP.py']:
 (O/p.name).write_bytes(p.read_bytes())
js({'status':'completed','models_retrained':False,'new_candidates_selected':False,'source_model_untouched':True},'extra_validation.json')
print('Extra audit completed')
