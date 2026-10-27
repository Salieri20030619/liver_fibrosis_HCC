import pathlib,gzip,csv,json
import pandas as pd
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');out=[]
for slide,gsm in [(1,'GSM8855042'),(2,'GSM8855043')]:
 pre=f'{gsm}_Slide{slide}_';cols=next(csv.reader(gzip.open(r/(pre+'exprMat_file.csv.gz'),'rt')))[2:];ctrl=[x for x in cols if x.startswith('SystemControl')];neg=[x for x in cols if 'NegPrb' in x or 'Negative' in x];genes=[x for x in cols if x not in ctrl+neg]
 m=pd.read_csv(r/(pre+'metadata_file.csv.gz'),usecols=['fov','CenterX_global_px','CenterY_global_px']);f=pd.read_csv(r/(pre+'fov_positions_file.csv.gz'))
 v={'slide_source':gsm,'assayed_features':len(cols),'system_controls':len(ctrl),'negative_probes':len(neg),'remaining_gene_probes':len(genes),'gene_probe_names_with_slash':[x for x in genes if '/' in x],'metadata_fov_missing_coordinates':sorted(set(m.fov)-set(f.FOV)),'coordinate_fovs_without_cells':sorted(set(f.FOV)-set(m.fov)),'duplicate_FOV_coordinates':int(f.FOV.duplicated().sum()),'missing_global_cell_coordinates':int(m[['CenterX_global_px','CenterY_global_px']].isna().any(axis=1).sum()),'template_gene_in_panel':{x:x in cols for x in ['DTNA','RUNX2','PLG','PARD3']}}
 out.append(v)
(r/'GSE292268_panel_coordinate_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
