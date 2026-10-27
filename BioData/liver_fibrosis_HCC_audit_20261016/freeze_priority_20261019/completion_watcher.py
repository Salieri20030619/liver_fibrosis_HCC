from pathlib import Path
import json,time,os,tarfile,csv,collections,datetime,traceback
R=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/freeze_priority_20261019')
def save(name,d):(R/name).write_text(json.dumps(d,ensure_ascii=False,indent=2))
def check(acc):
 state=json.loads((R/(acc+'_state.json')).read_text());rows=json.loads((R/(acc+'_matrix_rows.json')).read_text());clinical=json.loads((R/(acc+'_clinical_map.json')).read_text());by={x['GSM']:x for x in clinical};official=json.loads((R/(acc+'_official_filelist.json')).read_text());expected={x['Name']:int(x['Size']) for x in official if x['#Archive/File']=='File'}
 with tarfile.open(R/(acc+'_RAW.tar')) as t:actual={x.name:x.size for x in t.getmembers() if x.isfile()}
 members_ok=actual==expected
 joined=[]
 for x in rows:
  c=by.get(x['gsm']);y={**x,'clinical_join':c is not None}
  if c:
   y['donor']=c['donor'] if acc=='GSE289173' else c['patient_id'];y['disease_group']=c['disease group'] if acc=='GSE289173' else c['disease status'];y['clinical_conflict']=c.get('conflict_fields','');y['clinical_fields']=c
  joined.append(y)
 save(acc+'_matrix_clinical_join.json',joined)
 raw=[x for x in joined if acc=='GSE289173' or x.get('kind')=='raw'];donors={x['donor']:x for x in raw if x['clinical_join']};numok=all(x.get('nonfinite',0)==0 and x.get('negative',0)==0 and x.get('noninteger',0)==0 and x.get('bad_rows',0)==0 and x.get('dimensions_match',True) and x['unique_genes'] and x['unique_barcodes'] for x in raw)
 result={'accession':acc,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'official_archive_bytes_match':state['official_size_match'],'archive_sha256':state['sha256'],'official_hash_available':False,'archive_members_match_official_name_and_size':members_ok,'archive_missing_members':sorted(set(expected)-set(actual)),'archive_extra_members':sorted(set(actual)-set(expected)),'matrix_samples':len(raw),'raw_cells_total':sum(x['cells'] for x in raw),'gene_counts':sorted(set(x['genes'] for x in raw)),'raw_numerical_and_ID_checks_pass':numok,'matrix_GSMs_match_clinical':set(x['gsm'] for x in raw)==set(by),'independent_source_donor_labels':len(donors),'donor_groups':dict(collections.Counter(x['disease_group'] for x in donors.values())),'clinical_conflicts':[{'GSM':x['gsm'],'donor':x['donor'],'fields':x['clinical_conflict']} for x in raw if x.get('clinical_conflict')],'limitations':['This is data integrity checking only; no biological QC filtering, integration or training.','GSE289173 donor-level exact F and cell-type annotations not present in supplied raw matrices; group labels retained, aggregate F counts not assigned to people.','GSE202379 P6 diagnosis/SAF/age conflict held out; end-stage diagnosis not silently converted to F4.','No official checksum was supplied in downloaded filelist; computed SHA256 is a local content fingerprint.','Zero-library prevalence and author-filtered-cell correspondence have not yet been checked.']}
 result['file_integrity_pass']=members_ok and numok and result['matrix_GSMs_match_clinical'] and state['official_size_match'];save(acc+'_final_integrity.json',result)
 return result
save('completion_watcher_state.json',{'status':'waiting_for_download_and_matrix_validation','pid':os.getpid(),'primary_worker_pid':766772,'scope':'only the two approved archives'})
done=set()
try:
 while True:
  for acc in ['GSE289173','GSE202379']:
   if acc in done:continue
   p=R/(acc+'_state.json')
   try:d=json.loads(p.read_text())
   except (FileNotFoundError,json.JSONDecodeError):continue
   if d.get('status')=='complete':check(acc);done.add(acc)
  if len(done)==2:save('completion_watcher_state.json',{'status':'complete','resources':sorted(done),'pid':os.getpid()});break
  if not Path('/proc/766772').exists():save('completion_watcher_state.json',{'status':'blocked_primary_worker_exited_before_both_complete','completed_resources':sorted(done),'pid':os.getpid()});break
  time.sleep(30)
except Exception:
 save('completion_watcher_state.json',{'status':'error','error':traceback.format_exc(),'pid':os.getpid()});raise
