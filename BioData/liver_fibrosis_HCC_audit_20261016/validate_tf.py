import pathlib,gzip,csv,json,datetime,hashlib,collections
import numpy as np
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');p=r/'GSE227265_TFMotifActivity.csv.gz.complete';assert p.exists() and (r/'GSE227265_TF_resume_v2/complete.json').exists()
with gzip.open(p,'rt') as f:
 rd=csv.reader(f);head=next(rd);cells=head[1:];ids=[];nonfinite=0;badwidth=0;negative=0;low=float('inf');high=float('-inf')
 for row in rd:
  ids.append(row[0]);badwidth+=len(row)!=len(head);a=np.array(row[1:],dtype=float);nonfinite+=int((~np.isfinite(a)).sum());negative+=int((a<0).sum());low=min(low,float(a.min()));high=max(high,float(a.max()))
v={'resource':'GSE227265','file':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'rows_TF_motifs':len(ids),'columns_cells':len(cells),'duplicate_motif_IDs':len(ids)-len(set(ids)),'duplicate_cell_IDs':len(cells)-len(set(cells)),'bad_width_rows':badwidth,'nonfinite_values':nonfinite,'negative_values_expected_for_activity':negative,'range':[low,high],'cell_barcode_suffix_counts':dict(collections.Counter(x.rsplit('-',1)[-1] for x in cells)),'structural_pass':not(badwidth or nonfinite or len(cells)!=len(set(cells))),'semantic_alignment_pass':False,'limitations':['Cell barcode numeric suffix is not proven donor identifier','No matching peak-matrix column-name file; equal dimensions alone cannot prove column identity or order','Activity scores are processed ChromVAR output, not integer counts'],'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(r/'GSE227265_TF_full_validation.json').write_text(json.dumps(v,indent=2));(r/'GSE227265_TF_cell_ID_order.tsv').write_text('column_1based\tcell_id\n'+'\n'.join(str(i)+'\t'+x for i,x in enumerate(cells,2)));print(json.dumps(v))
