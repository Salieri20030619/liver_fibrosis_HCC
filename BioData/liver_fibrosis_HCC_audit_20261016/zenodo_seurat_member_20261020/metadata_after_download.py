from pathlib import Path
import json,time,subprocess,datetime,hashlib
P=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/zenodo_seurat_member_20261020')
def save(v):
 v['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();q=P/'metadata_workflow_state.json.tmp';q.write_text(json.dumps(v,indent=2));q.replace(P/'metadata_workflow_state.json')
save({'status':'waiting_for_verified_member'})
while True:
 s=json.loads((P/'state.json').read_text())
 if s['status']=='blocked_or_failed':save({'status':'download_failed_metadata_not_started','download_error':s.get('error')});raise SystemExit(2)
 if s['status']=='download_and_member_integrity_complete':break
 time.sleep(20)
i=json.loads((P/'member_integrity.json').read_text())
assert i['bytes']==2442015643 and i['crc32']=='11738af9' and i['crc32_matches_zip_directory']
save({'status':'reading_metadata','source_sha256':i['sha256']})
with (P/'metadata_R_console.log').open('wb') as log:
 r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/r_seurat/bin/Rscript','--vanilla',str(P/'metadata_audit.R'),str(P)],stdout=log,stderr=subprocess.STDOUT)
if r.returncode:save({'status':'metadata_read_failed','returncode':r.returncode});raise SystemExit(r.returncode)
files=[]
for f in (P/'metadata_audit').iterdir():
 if not f.is_file():continue
 h=hashlib.sha256()
 with f.open('rb') as inp:
  for chunk in iter(lambda:inp.read(1024*1024),b''):h.update(chunk)
 files.append({'name':f.name,'bytes':f.stat().st_size,'sha256':h.hexdigest()})
save({'status':'metadata_extraction_complete','source_sha256':i['sha256'],'files':files})
