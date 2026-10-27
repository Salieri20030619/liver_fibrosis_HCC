import pathlib,json,datetime,hashlib
import pandas as pd,numpy as np
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');o=r/'CosMx_Slide1_QC_coordinates';o.mkdir(exist_ok=False)
rp=r/'GSM8855042_RDS_inspection/cell_metadata.tsv';cp=r/'GSM8855042_Slide1_metadata_file.csv.gz'
a=pd.read_csv(rp,sep='\t');b=pd.read_csv(cp).set_index('cell').loc[a.cell_id].reset_index();assert a.cell_id.tolist()==b.cell.tolist()
qc=[c for c in a if c.startswith('qcFlags')];v={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cells':len(a),'QC':{},'coordinates':{},'scope':'Existing supplied QC flags only; no filtering, clustering or patient assignment. Geometric affine coefficients below are empirical numeric equivalence checks, not biological modelling.'}
for c in qc:v['QC'][c]={'counts':a[c].fillna('<missing>').value_counts().to_dict(),'CSV_flag_mismatches':int((a[c].fillna('<missing>')!=b[c].fillna('<missing>')).sum())}
joint=a[qc].eq('Pass').all(axis=1);v['joint_all_six_pass']=int(joint.sum());v['joint_not_all_pass']=int((~joint).sum());v['joint_cell_flags_pass']=int(a[[c for c in qc if c!='qcFlagsFOV']].eq('Pass').all(axis=1).sum())
v['FOV_flag_within_FOV_consistency']=int((a.groupby('fov').qcFlagsFOV.nunique()!=1).sum());v['FOV_level_counts']=a[['fov','qcFlagsFOV']].drop_duplicates().qcFlagsFOV.value_counts().to_dict()
lab=next(c for c in a if c.endswith('_clusters'));un=a[lab].isin(list('abcdefghij'));v['letter_clusters_cells']=int(un.sum());v['letter_clusters_percent']=float(un.mean()*100);v['letter_clusters_joint_QC_pass']=int((un&joint).sum());v['named_labels_joint_QC_pass']=int((~un&joint).sum())
for rc,cc in [('x_FOV_px','CenterX_local_px'),('y_FOV_px','CenterY_local_px')]:
 delta=a[rc].to_numpy()-b[cc].to_numpy();v['coordinates'][rc]={'CSV_column':cc,'max_absolute_error_px':float(np.max(np.abs(delta))),'mismatches':int(np.count_nonzero(delta))}
for rc,cc in [('x_slide_mm','CenterX_global_px'),('y_slide_mm','CenterY_global_px')]:
 x=b[cc].to_numpy(float);y=a[rc].to_numpy(float);mask=np.isfinite(x)&np.isfinite(y);slope,offset=np.polyfit(x[mask],y[mask],1);res=y-(slope*x+offset)
 scale=0.000120280945;direct=y-scale*x
 v['coordinates'][rc]={'CSV_column':cc,'finite_pairs':int(mask.sum()),'empirical_mm_per_pixel':float(slope),'empirical_offset_mm':float(offset),'affine_max_abs_error_mm':float(np.max(np.abs(res[mask]))),'affine_error_gt_1e-9_mm':int((np.abs(res[mask])>1e-9).sum()),'tested_fixed_mm_per_pixel':scale,'fixed_scale_no_offset_max_abs_error_mm':float(np.max(np.abs(direct[mask]))),'fixed_scale_no_offset_error_gt_1e-9_mm':int((np.abs(direct[mask])>1e-9).sum())}
fv=pd.read_csv(r/'CosMx_Slide1_mapping_evidence/latest_fovs.tsv',sep='\t');fc=pd.read_csv(r/'GSM8855042_Slide1_fov_positions_file.csv.gz');j=fv.merge(fc,on=['Slide','FOV'],suffixes=('_rds','_csv'),validate='one_to_one');v['FOV_tables']={'rds_shape':list(fv.shape),'csv_shape':list(fc.shape),'matched_rows':len(j),'coordinate_max_abs_differences':{c:float(np.max(np.abs(j[c+'_rds']-j[c+'_csv']))) for c in ['X_mm','Y_mm']},'patient_mapping':'Still absent; coordinates not used to infer patients'}
pd.DataFrame([{'flag':c,**z['counts'],'CSV_mismatches':z['CSV_flag_mismatches']} for c,z in v['QC'].items()]).to_csv(o/'QC_flag_counts.tsv',sep='\t',index=False)
a.groupby(qc,dropna=False).size().reset_index(name='cells').to_csv(o/'QC_joint_patterns.tsv',sep='\t',index=False)
v['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [rp,cp]}
with (o/'QC_coordinate_audit.json').open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
print(json.dumps(v))
