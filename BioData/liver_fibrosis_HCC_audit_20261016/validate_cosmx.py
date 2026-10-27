import pathlib,json,gzip,hashlib,datetime
import pandas as pd,numpy as np
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');out=[]
for slide,gsm in [(1,'GSM8855042'),(2,'GSM8855043')]:
 pre=f'{gsm}_Slide{slide}_';files={k:r/(pre+v+'.csv.gz') for k,v in [('expression','exprMat_file'),('metadata','metadata_file'),('fov','fov_positions_file')]}
 if not all(p.exists() for p in files.values()):continue
 e=pd.read_csv(files['expression']);m=pd.read_csv(files['metadata']);f=pd.read_csv(files['fov'])
 ek=e[['fov','cell_ID']].astype(str).agg('_'.join,axis=1);mk=m[['fov','cell_ID']].astype(str).agg('_'.join,axis=1);a=e.drop(columns=['fov','cell_ID']).to_numpy()
 v={'slide':slide,'GSM':gsm,'expression_cells':len(e),'metadata_cells':len(m),'assayed_features':e.shape[1]-2,'duplicate_expression_keys':int(ek.duplicated().sum()),'duplicate_metadata_keys':int(mk.duplicated().sum()),'expression_only_cells':len(set(ek)-set(mk)),'metadata_only_cells':len(set(mk)-set(ek)),'same_cell_order':ek.tolist()==mk.tolist(),'nonfinite':int((~np.isfinite(a)).sum()),'negative':int((a<0).sum()),'noninteger':int((a!=np.floor(a)).sum()),'fov_rows':len(f),'fov_columns':f.columns.tolist(),'metadata_fovs':sorted(m.fov.unique().tolist()),'fov_data':f.to_dict('records'),'metadata_columns':m.columns.tolist(),'sample_fields':{c:m[c].value_counts(dropna=False).to_dict() for c in m if any(x in c.lower() for x in ['sample','patient','donor','tissue','slide'])},'files':{k:{'file':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in files.items()}}
 pd.DataFrame({'slide':slide,'key_fov_cellID':ek,'in_metadata':ek.isin(mk)}).to_csv(r/f'GSE292268_Slide{slide}_cell_alignment.tsv',sep='\t',index=False)
 out.append(v)
(r/'GSE292268_full_validation.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'slides':out},indent=2,default=str))
print(json.dumps([{k:v for k,v in x.items() if k not in ['metadata_columns','fov_data','files','metadata_fovs']} for x in out],default=str))
