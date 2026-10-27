import pathlib,json,time,os,signal,subprocess,datetime,importlib.util,csv,gzip,tarfile,hashlib,collections
R=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/round2_processed_20261017')
def save(name,obj):(R/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
small=['GSE246221_rawcounts_allsamples_human.txt.gz','GSE246221_rawcounts_allsamples_mouse.txt.gz','GSE83596_series_matrix.txt.gz'];archive='GSE244832_hLIVER_processed_files.tar.gz'
save('range_controller_state.json',{'pid':os.getpid(),'status':'waiting_for_small_validations','started_utc':now()})
deadline=time.time()+3600
while time.time()<deadline:
 states=[json.loads((R/(s+'.state.json')).read_text()) for s in small]
 if all(v.get('status')=='completed' for v in states):break
 if any(v.get('status') in ['error','download_incomplete'] for v in states):raise RuntimeError('Small-file validation needs review; original job left running')
 time.sleep(30)
else:raise TimeoutError('Small downloads not finished; original job left running')
pid=json.loads((R/'job_identity.json').read_text())['pid'];proc=pathlib.Path('/proc')/str(pid)/'cmdline'
if proc.exists():
 assert b'round2_server_downloads.py' in proc.read_bytes()
 os.killpg(pid,signal.SIGTERM);time.sleep(2)
p=R/archive;prefix=R/(archive+'.part')
save('range_controller_state.json',{'pid':os.getpid(),'status':'range_resume_running','prefix_bytes':prefix.stat().st_size if prefix.exists() else 0,'utc':now(),'small_files_completed':True,'concurrency':12})
if not p.exists():
 cmd=['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python',str(R/'resume_verified_direct.py'),'GSE244832',archive,'727189398',str(R),str(prefix),'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244832/suppl/'+archive]
 with (R/'archive_range_resume.log').open('ab') as log:cp=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 completed=R/(archive+'.complete')
 if cp.returncode!=0 or not completed.exists():
  save('range_controller_state.json',{'pid':os.getpid(),'status':'range_download_incomplete','returncode':cp.returncode,'utc':now()});raise RuntimeError('Range download incomplete; verified pieces preserved')
 completed.replace(p)
spec=importlib.util.spec_from_file_location('audit_helpers',R/'round2_server_downloads.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
state=json.loads((R/(archive+'.state.json')).read_text());state.update(status='download_complete_validating',bytes=p.stat().st_size,sha256=mod.digest(p),range_resume=True);save(archive+'.state.json',state)
state['validation']=mod.validate(p,'archive');state['status']='completed';state['finished_utc']=now();save(archive+'.state.json',state)
# Retain all small members; inventory tabular fields without assuming assay pairing.
folder=R/'complete_small_members';folder.mkdir(exist_ok=True);inventory=[]
with tarfile.open(p,'r:gz') as tf:
 for m in tf:
  if not m.isfile() or m.size>100000000:continue
  dest=folder/(hashlib.sha256(m.name.encode()).hexdigest()[:8]+'_'+pathlib.PurePosixPath(m.name).name)
  if not dest.exists():
   with tf.extractfile(m) as f,dest.open('xb') as w:
    for b in iter(lambda:f.read(1048576),b''):w.write(b)
  rec={'member':m.name,'bytes':m.size,'path':str(dest),'sha256':mod.digest(dest)}
  if any(x in m.name.lower() for x in ['csv','tsv','meta','gene','cell']):
   try:
    magic=dest.open('rb').read(2);op=gzip.open if magic==b'\x1f\x8b' else open
    with op(dest,'rt',encoding='utf-8',newline='') as f:
     first=f.readline();delimiter='\t' if '\t' in first else ',';head=next(csv.reader([first],delimiter=delimiter));reader=csv.reader(f,delimiter=delimiter)
     n=0;sets=[set() for _ in head];examples=[]
     for row in reader:
      n+=1
      if len(examples)<3:examples.append(row)
      for i,x in enumerate(row[:len(head)]):
       if len(sets[i])<=200:sets[i].add(x)
     rec.update(columns=head,rows=n,examples=examples,field_values={head[i]:sorted(s) if len(s)<=200 else 'more than 200 distinct values' for i,s in enumerate(sets)})
   except Exception as e:rec['parse_error']=repr(e)
  inventory.append(rec)
save('complete_small_member_field_audit.json',inventory)
save('range_controller_state.json',{'pid':os.getpid(),'status':'completed','utc':now(),'archive_validated':True,'small_member_files':len(inventory)})
print(json.dumps({'status':'completed','archive_bytes':p.stat().st_size,'metadata_files':len(inventory)}),flush=True)
