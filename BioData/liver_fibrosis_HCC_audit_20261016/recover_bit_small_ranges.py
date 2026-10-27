import pathlib,json,subprocess,datetime,time,hashlib,os,fcntl,shutil
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016')
lock=(r/'BIT_range_resume/recovery.lock').open('a')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
for line in subprocess.check_output(['ps','-eo','pid,args'],text=True).splitlines():
 if str(r)+'/bit_resumable_extract.py' in line or str(r)+'/finish_download_audit.py' in line:raise RuntimeError('Original BIT process still active: '+line)
source=(r/'bit_resumable_extract.py').read_text()
replacement='''def fetch(i):
 a=start+i*chunk;b=min(total-1,a+chunk-1);log=d/f'{i:05}.json';part=d/f'{i:05}.part'
 if log.exists() and part.exists():
  x=json.loads(log.read_text())
  if x.get('verified') and x['range']==[a,b] and part.stat().st_size==b-a+1:return part
  raise RuntimeError('Existing canonical part inconsistent; preserve for review')
 assert not log.exists() and not part.exists()
 pieces=[]
 for j,aa in enumerate(range(a,b+1,1048576)):
  bb=min(b,aa+1048575);base=d/f'{i:05}.sub{j:02}';sub=pathlib.Path(str(base)+'.part');meta=pathlib.Path(str(base)+'.json')
  if sub.exists() and meta.exists():
   v=json.loads(meta.read_text());assert v['verified'] and v['range']==[aa,bb] and sub.stat().st_size==bb-aa+1
   pieces.append(sub);continue
  assert not sub.exists() and not meta.exists()
  for attempt in range(3):
   if attempt:time.sleep(60*2**(attempt-1))
   endpoint=url if attempt<2 else 'https://zenodo.org/records/14231098/files/Input_Data.tar.gz?download=1'
   tag=str(time.time_ns());h=pathlib.Path(str(base)+'.'+tag+'.headers');tmp=pathlib.Path(str(base)+'.'+tag+'.partial');beg=datetime.datetime.now(datetime.timezone.utc).isoformat()
   p=subprocess.run(['curl','-sS','-L','--fail','--connect-timeout','30','--max-time','300','--max-filesize',str(bb-aa+1),'--range',f'{aa}-{bb}','-D',str(h),'-o',str(tmp),endpoint],capture_output=True,text=True)
   hs=h.read_text() if h.exists() else '';n=tmp.stat().st_size if tmp.exists() else 0
   ok=p.returncode==0 and n==bb-aa+1 and f'bytes {aa}-{bb}/{total}' in hs
   v={'url':endpoint,'range':[aa,bb],'bytes':n,'verified':ok,'curl_exit':p.returncode,'stderr':p.stderr,'http_status_lines':[s for s in hs.splitlines() if s.startswith('HTTP/')],'headers':str(h),'started_utc':beg,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
   with pathlib.Path(str(base)+'.'+tag+'.json').open('x') as f:json.dump(v,f,indent=2)
   if ok:
    v['sha256']=hashlib.sha256(tmp.read_bytes()).hexdigest();os.link(tmp,sub)
    with meta.open('x') as f:json.dump(v,f,indent=2)
    pieces.append(sub);print(json.dumps({'part':i,'subpart':j,'verified_new_bytes':n}),flush=True);break
  else:raise RuntimeError(f'Bounded smaller-range retries failed {aa}-{bb}; retained all valid parts')
 with part.open('xb') as f:
  for sub in pieces:
   with sub.open('rb') as src:shutil.copyfileobj(src,f)
 assert part.stat().st_size==b-a+1
 v={'url':url,'range':[a,b],'bytes':part.stat().st_size,'verified':True,'sha256':hashlib.sha256(part.read_bytes()).hexdigest(),'method':'Four or fewer exact HTTP 206 subranges, 1 MiB each','subpart_metadata':[str(p.with_suffix('.json')) for p in pieces],'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 with log.open('x') as f:json.dump(v,f,indent=2)
 print(json.dumps({'part':i,'bytes':v['bytes'],'verified':True}),flush=True);return part
'''
source=source[:source.index('def fetch(i):')]+replacement+source[source.index('pool=concurrent.futures.ThreadPoolExecutor'):]
exec(compile(source,'bit_small_range_recovery','exec'),globals())
if (r/'BIT_range_resume/Liver.h5ad.complete').exists() and not (r/'BIT_h5_inspection').exists():
 subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python',str(r/'inspect_bit_h5.py')],check=True)
