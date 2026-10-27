from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad,scanpy as sc,traceback
P=Path(__file__).parent;O=P/'continuation_doublet_diagnostics_v1';O.mkdir(exist_ok=False);a=ad.read_h5ad(P/'continuation_annotation_evidence_v1/core55_record_counts.h5ad',backed='r');m=a.obs.copy();m['source_record']=m.study.astype(str)+'|'+m.donor.astype(str);rows=[];tables=[]
for record in sorted(m.source_record.unique()):
 try:
  s=a[m.source_record.eq(record).to_numpy()].to_memory();sc.pp.scrublet(s,expected_doublet_rate=.06,n_prin_comps=30,random_state=20261020,verbose=False);d=s.obs[['doublet_score','predicted_doublet']].copy();d['source_record']=record;tables.append(d);rows.append({'source_record':record,'cells':len(s),'called_doublets':int(d.predicted_doublet.sum()),'threshold':float(s.uns['scrublet']['threshold']),'status':'completed'})
 except Exception as e:rows.append({'source_record':record,'status':'failed','error':repr(e)});(O/(record.replace('|','__')+'_error.txt')).write_text(traceback.format_exc())
 pd.DataFrame(rows).to_csv(O/'record_diagnostics.tsv',sep='\t',index=False);print(rows[-1],flush=True)
pd.concat(tables).to_csv(O/'cell_doublet_diagnostics.tsv.gz',sep='\t');(O/'validation.json').write_text(json.dumps({'status':'completed' if len(tables)==55 else 'partial','records_attempted':55,'records_completed':len(tables),'expected_doublet_rate':.06,'parameter_is_reconstructed':True,'interpretation':'Per-source-record Scrublet simulation diagnostic; nuclei/scRNA source heterogeneity and unknown loading limits calibration. No automatic removal or claim that marker co-detection proves doublets. No empty droplets for ambient estimation.'},indent=2))
