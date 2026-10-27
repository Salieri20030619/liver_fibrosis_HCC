from pathlib import Path
import json,hashlib,datetime,tarfile,csv
P=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/end_to_end_20261020');O=P/'ml113_repaired_20261020_v1'
assert (O/'results/validation.json').exists()
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');B=P/('ml113_delivery_'+stamp);B.mkdir()
selected=[];rows=[]
for d in ['code','inputs','results','code_revision_api_v2','code_revision_compat_v3','code_revision_PLS_v4','jobs']:
 for q in (O/d).rglob('*'):
  if not q.is_file():continue
  size=q.stat().st_size
  include=(d!='jobs' or q.suffix in ['.tsv','.json','.txt'] or q.name in ['preprocessing.rds','selected_frozen_model.rds']) and size<25000000
  h=hashlib.sha256(q.read_bytes()).hexdigest() if size<50000000 else ''
  rows.append({'path':str(q.relative_to(O)),'bytes':size,'sha256':h,'included':include})
  if include:selected.append(q)
for q in O.iterdir():
 if q.is_file() and q.suffix in ['.json','.tsv','.txt','.R','.log']:
  rows.append({'path':q.name,'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'included':True});selected.append(q)
for q in (O/'package_sources').iterdir():
 if q.is_file():rows.append({'path':str(q.relative_to(O)),'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'included':False})
with open(B/'manifest.tsv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes','sha256','included'],delimiter='\t');w.writeheader();w.writerows(rows)
arc=B/'ml113_audit_code_results.tar.gz'
with tarfile.open(arc,'w:gz') as t:
 for q in selected:t.add(q,arcname=q.relative_to(O))
 t.add(B/'manifest.tsv',arcname='artifact_manifest.tsv')
 for q in P.glob('cont_ml113*.py'):t.add(q,arcname='execution_helpers/'+q.name)
v={'utc':stamp,'archive':str(arc),'bytes':arc.stat().st_size,'sha256':hashlib.sha256(arc.read_bytes()).hexdigest(),'included_files':len(selected),'source_root':str(O),'limits':'Large inner prediction RDS and installed package binaries remain on server. File manifests and all113 configuration metrics included; no formalfigures.'}
(B/'receipt.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
