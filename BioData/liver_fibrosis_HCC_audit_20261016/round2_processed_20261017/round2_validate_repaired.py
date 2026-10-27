import pathlib,json,time,os,subprocess,datetime,importlib.util,csv,gzip,tarfile,hashlib,collections
R=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/round2_processed_20261017')
archive='GSE244832_hLIVER_processed_files.tar.gz';p=R/archive
def save(name,obj):(R/name).write_text(json.dumps(obj,indent=2))
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
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
