import os,json,time,signal,threading,subprocess,hashlib,re,shutil,concurrent.futures,datetime,fcntl
from pathlib import Path
from verified_chunk_io import append_verified
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
lock=open(ROOT/'download.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
PLAN=json.loads((ROOT/'manifest.json').read_text())
assert PLAN[0]['filename']=='NH4_S2_L004_R2_001.fastq.gz'
assert PLAN[0]['url']=='https://ftp.sra.ebi.ac.uk/vol1/run/ERR116/ERR11656833/NH4_S2_L004_R2_001.fastq.gz'
assert PLAN[0]['md5']=='c758d29c0cacc9b720eac8cd95a2d507'
assert PLAN[0]['initial_bytes']==0
assert Path(PLAN[0]['destination']).resolve()==(ROOT/'data').resolve()
assert not (ROOT/'data'/PLAN[0]['filename']).exists()
assert not (ROOT/'data'/(PLAN[0]['filename']+'.partial')).exists()
prior=ROOT.parent/'hammond_connection_retry_20261025_1400'
assert json.loads((prior/'progress.json').read_text())['finished']
assert not Path('/proc/'+str(json.loads((prior/'launch.json').read_text())['pid'])).exists()
prior_lock=open(prior/'download.lock','a');fcntl.flock(prior_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)

assert len(PLAN)==1 and sum(x['remaining_bytes'] for x in PLAN)==3900472193
assert len({x['filename'] for x in PLAN})==1
e_lock=open(ROOT.parent/'hammond_e_download_20261023'/'download.lock','a');fcntl.flock(e_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
STOP=ROOT.parent/'hammond_overnight_20261024'/'STOP'; halt=threading.Event(); mutex=threading.Lock();states={};children={}
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(path,obj):
 tmp=Path(str(path)+'.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2));os.replace(tmp,path)
def event(obj):
 with mutex:
  with (ROOT/'events.jsonl').open('a') as f:f.write(json.dumps({'utc':utc(),**obj})+'\n')
def stopped():return STOP.exists() or (ROOT.parent/'hammond_e_download_20261023'/'STOP').exists() or halt.is_set()
def sigstop(sig,frame):
 STOP.touch();halt.set()
signal.signal(signal.SIGTERM,sigstop);signal.signal(signal.SIGINT,sigstop)
if STOP.exists():raise SystemExit('STOP exists: no automatic restart permitted')
if shutil.disk_usage(ROOT).free<sum(int(x['bytes']) for x in PLAN)+1024**3:raise SystemExit('Insufficient conservative disk reserve')
write(ROOT/'process.json',{'pid':os.getpid(),'started':utc(),'proxy':'http://127.0.0.1:17897','max_concurrent':1,'total_bytes':sum(int(x['bytes']) for x in PLAN),'initial_bytes':sum(x['initial_bytes'] for x in PLAN),'stop_file':str(STOP),'auto_restart':False})
def md5(path):
 h=hashlib.md5()
 with path.open('rb') as f:
  while True:
   if stopped():raise InterruptedError('Stopped during MD5')
   b=f.read(8*1024*1024)
   if not b:break
   h.update(b)
 return h.hexdigest()
def state(name,**kw):
 with mutex:states.setdefault(name,{}).update(kw)
def worker(item):
 name=item['filename'];size=int(item['bytes']);dest=Path(item['destination']);dest.mkdir(parents=True,exist_ok=True);final=dest/name;part=dest/(name+'.partial')
 state(name,state='checking_existing',expected_bytes=size,committed_bytes=item['initial_bytes'],active_chunk_bytes=0)
 try:
  if stopped():raise InterruptedError('User stop before file start')
  if final.exists():
   if final.stat().st_size==size and md5(final)==item['md5']:
    state(name,state='complete_verified',committed_bytes=size,md5=item['md5']);return
   raise RuntimeError('Existing final file fails size/MD5; preserved without overwrite')
  offset=part.stat().st_size if part.exists() else 0
  if offset!=item['initial_bytes']:raise RuntimeError('Partial changed since preflight; refuse uncertain resume')
  # All committed bytes are exact validated contiguous ranges; no test fragments imported.
  while offset<size:
   if stopped():raise InterruptedError('User stop')
   end=min(offset+1*1024*1024,size)-1;want=end-offset+1;ok=False
   for attempt in range(1,4):
    if stopped():raise InterruptedError('Stop requested')
    stem=ROOT/'chunks'/f'{name}.{offset}.{attempt}'
    body=Path(str(stem)+'.partial');headers=Path(str(stem)+'.headers');stderr=Path(str(stem)+'.stderr')
    if body.exists():raise RuntimeError('Existing chunk attempt preserved; manual inspection required')
    cmd=['curl','--disable','--proxy','http://127.0.0.1:17897','--noproxy','','--silent','--show-error','--fail','--connect-timeout','20','--max-time','180','--max-filesize',str(want),'--range',f'{offset}-{end}','--dump-header',str(headers),'--output',str(body),item['url']]
    state(name,state='downloading',committed_bytes=offset,active_chunk=str(body),range_start=offset,range_end=end,attempt=attempt)
    with stderr.open('wb') as ef:
     p=subprocess.Popen(cmd,stderr=ef,stdout=subprocess.DEVNULL,start_new_session=True)
     with mutex:children[name]=p.pid
     while p.poll() is None:
      state(name,active_chunk_bytes=body.stat().st_size if body.exists() else 0)
      if stopped():
       p.terminate()
       try:p.wait(timeout=3)
       except subprocess.TimeoutExpired:p.kill()
       p.wait();raise InterruptedError('Stopped; active fragment retained')
      time.sleep(.5)
    with mutex:children.pop(name,None)
    h=headers.read_text(errors='replace') if headers.exists() else ''
    codes=re.findall(r'^HTTP/\S+\s+(\d+)',h,re.M);code=int(codes[-1]) if codes else 0
    ranges=re.findall(r'^content-range:\s*(.+)$',h,re.M|re.I);cr=ranges[-1].strip() if ranges else None
    lengths=re.findall(r'^content-length:\s*(\d+)',h,re.M|re.I)
    actual=body.stat().st_size if body.exists() else 0
    valid=p.returncode==0 and code==206 and cr==f'bytes {offset}-{end}/{size}' and actual==want and (not lengths or int(lengths[-1])==want)
    event({'file':name,'range':[offset,end],'attempt':attempt,'http':code,'content_range':cr,'received_bytes':actual,'curl_rc':p.returncode,'accepted':valid})
    if valid:
     body_hash=hashlib.sha256(body.read_bytes()).hexdigest()
     event({'file':name,'action':'before_append','offset':offset,'bytes':want,'body_sha256':body_hash})
     proof=append_verified(body,part,offset,want)
     if proof['body_sha256']!=body_hash:raise RuntimeError('Body changed during append; preserve evidence')
     event({'file':name,'action':'persisted_chunk_verified',**proof})
     body.unlink();offset+=want;state(name,committed_bytes=offset,active_chunk_bytes=0);ok=True;break
    if code in [401,403,407,429]:raise RuntimeError(f'HTTP {code}; no bypass or retry')
    if p.returncode not in [5,6,7,18,28,35,52,55,56] and code not in [500,502,503,504]:raise RuntimeError('Nontransient response/range/length failure')
    if attempt<3:
     for _ in range(40 if attempt==1 else 120):
      if stopped():raise InterruptedError('Stop during bounded retry delay')
      time.sleep(.5)
   if not ok:raise RuntimeError('Transient error retry limit reached (3 total attempts)')
  state(name,state='verifying_MD5',committed_bytes=offset,active_chunk_bytes=0)
  digest=md5(part)
  event({'file':name,'action':'full_md5','expected':item['md5'],'observed':digest,'bytes':part.stat().st_size})
  if digest!=item['md5']:raise RuntimeError('Full-file MD5 mismatch; partial preserved')
  if stopped():raise InterruptedError('Stop before final rename')
  os.replace(part,final);state(name,state='complete_verified',md5=digest,committed_bytes=size)
  event({'file':name,'complete_verified':True,'bytes':size,'md5':digest})
 except InterruptedError as e:state(name,state='stopped',error=str(e),active_chunk_bytes=0)
 except Exception as e:
  state(name,state='failed',error=repr(e),active_chunk_bytes=0);event({'file':name,'error':repr(e),'continue_other_files':True})
 finally:
  with mutex:children.pop(name,None)
def snapshot(done=False):
 with mutex:d={'utc':utc(),'pid':os.getpid(),'finished':done,'stop_flag':STOP.exists(),'halt':halt.is_set(),'files':json.loads(json.dumps(states)),'child_pids':dict(children)}
 d['verified_files']=sum(v.get('state')=='complete_verified' for v in d['files'].values());d['expected_files']=len(PLAN);d['expected_bytes']=sum(int(x['bytes']) for x in PLAN)
 d['initial_bytes']=sum(x['initial_bytes'] for x in PLAN);d['remaining_at_launch']=3900472193
 d['failed_files']=[k for k,v in d['files'].items() if v.get('state')=='failed']
 d['committed_bytes']=sum(v.get('committed_bytes',0) for v in d['files'].values());d['inflight_bytes']=sum(v.get('active_chunk_bytes',0) for v in d['files'].values())
 d['new_committed_bytes']=d['committed_bytes']-d['initial_bytes']
 write(ROOT/'progress.json',d)
(ROOT/'chunks').mkdir(exist_ok=True)
event({'action':'start','pid':os.getpid(),'scope':'independent MD5 repair c/NH4_S2_L004_R2_001.fastq.gz only'})
for x in PLAN:state(x['filename'],state='queued',expected_bytes=int(x['bytes']),committed_bytes=x['initial_bytes'],active_chunk_bytes=0)
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
 fs=[pool.submit(worker,x) for x in PLAN]
 while not all(f.done() for f in fs):snapshot();time.sleep(5)
 for f in fs:f.result()
snapshot(True);event({'action':'finished','user_stop':STOP.exists()})
