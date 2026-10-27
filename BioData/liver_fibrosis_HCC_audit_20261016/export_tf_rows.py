import pathlib,gzip,csv,json,hashlib
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');p=r/'GSE227265_TFMotifActivity.csv.gz.complete';ids=[]
with gzip.open(p,'rt') as f:
 rd=csv.reader(f);next(rd)
 for row in rd:ids.append(row[0])
out=r/'GSE227265_TF_motif_row_IDs.tsv'
with out.open('x') as f:f.write('row_1based\tmotif_id\n'+'\n'.join(str(i)+'\t'+x for i,x in enumerate(ids,2)))
print(json.dumps({'rows':len(ids),'examples':ids[:20],'RUNX_related_rows':[x for x in ids if 'RUNX' in x.upper()],'file':str(out),'bytes':out.stat().st_size,'sha256':hashlib.sha256(out.read_bytes()).hexdigest()}))
