import pathlib,json,datetime,hashlib,zipfile,csv
import pandas as pd
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');d=r/'CosMx_Slide2_mapping_evidence'
a=pd.read_csv(d/'cell_ID_fields.tsv',sep='\t',dtype=str);m=pd.read_csv(r/'GSM8855043_Slide2_metadata_file.csv.gz',dtype=str);e=pd.read_csv(r/'GSM8855043_Slide2_exprMat_file.csv.gz',usecols=['fov','cell_ID'],dtype=str)
ri=a.rds_row_id.tolist();ci=m['cell'].tolist();assert len(set(ci))==len(ci)
cm=m.set_index('cell');aligned=cm.loc[ri]
exprkeys=e[['fov','cell_ID']].agg('_'.join,axis=1);metakeys=m[['fov','cell_ID']].agg('_'.join,axis=1)
join=pd.DataFrame({'rds_row_1based':range(1,len(ri)+1),'cell_id':ri,'CSV_row_1based':pd.Series(range(1,len(ci)+1),index=ci).loc[ri].to_numpy(),'FOV':a.fov,'csv_fov':aligned.fov.to_numpy(),'csv_cell_ID':aligned.cell_ID.to_numpy()})
join.to_csv(d/'RDS_CSV_exact_cell_join.tsv',sep='\t',index=False)
f=pd.read_csv(d/'latest_fovs.tsv',sep='\t');ann=pd.read_csv(d/'selected_cell_annotations.tsv',sep='\t',dtype=str)
fields={c:ann[c].value_counts(dropna=False).head(10).to_dict() for c in ann.columns if any(s in c.lower() for s in ['patient','donor','sample','tissue','slide','orig.ident'])}
clusters=[c for c in ann if c.endswith('_clusters')];counts=[]
for c in clusters:
 for label,n in ann[c].value_counts(dropna=False).items():counts.append({'field':c,'label':label,'cells':int(n)})
pd.DataFrame(counts).to_csv(d/'provided_celltype_label_counts.tsv',sep='\t',index=False)
v={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'resource':'GSM8855043 / GSE292268 Slide2','RDS_provenance':json.loads((r/'GSM8855043_RDS_resume_v2/complete.json').read_text()),'rds_cells':len(ri),'CSV_cells':len(ci),'RDS_duplicate_IDs':len(ri)-len(set(ri)),'CSV_duplicate_IDs':len(ci)-len(set(ci)),'same_ID_set':set(ri)==set(ci),'same_ID_order':ri==ci,'RDS_cell_id_field_equals_rownames':a.cell_id.tolist()==ri,'RDS_cell_ID_field_equals_rownames':a.cell_ID.tolist()==ri,'FOV_mismatches_after_exact_ID_join':int((a.fov.to_numpy()!=aligned.fov.to_numpy()).sum()),'CSV_expression_metadata_FOV_cellID_same_order':exprkeys.tolist()==metakeys.tolist(),'assay_ID_alignment':pd.read_csv(d/'assay_ID_alignment.tsv',sep='\t').to_dict('records'),'FOV_rows':len(f),'FOV_unique':int(f.FOV.nunique()),'FOV_fields':{c:{'distinct':int(f[c].nunique()),'values':f[c].unique().tolist() if f[c].nunique()<10 else None} for c in f},'sample_like_metadata_fields':fields,'patient_mapping_status':'Not recovered: no explicit patient/donor/biopsy assignment in inspected meta.data, misc latest.fovs, tools or images. Constant Run_Tissue_name/slide identifiers cannot assign the two biopsies. Do not infer patient from FOV ordering or coordinates.','provided_annotations':'30 existing celltype/cluster labels include letter-labelled unresolved clusters; these are supplied annotations, not newly validated cell types. Existing niche labels and QC flags retained.','scope':'Object metadata and assay ID alignment only; expression values and biological annotation accuracy not newly revalidated.'}
with (d/'Slide2_annotation_audit.json').open('x') as out:json.dump(v,out,indent=2,ensure_ascii=False)
print(json.dumps(v,ensure_ascii=False))
