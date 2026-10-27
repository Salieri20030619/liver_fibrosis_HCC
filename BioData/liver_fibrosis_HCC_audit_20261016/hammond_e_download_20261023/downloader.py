import os,json,time,signal,threading,subprocess,hashlib,re,shutil,concurrent.futures,datetime,fcntl
from pathlib import Path
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
lock=open(ROOT/'download.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
PLAN=json.loads((ROOT/'manifest.json').read_text())
assert len(PLAN)==4 and {x['run'] for x in PLAN}=={'ERR14761304'}
assert sum(int(x['bytes']) for x in PLAN)==6729410219
STOP=ROOT/'STOP'; halt=threading.Event(); mutex=threading.Lock();states={};children={}
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(path,obj):
 tmp=Path(str(path)+'.tmp');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2));os.replace(tmp,path)
def event(obj):
 with mutex:
  with (ROOT/'events.jsonl').open('a') as f:f.write(json.dumps({'utc':utc(),**obj})+'\n')
def stopped():return STOP.exists() or halt.is_set()
def sigstop(sig,frame):
 STOP.touch();halt.set()
signal.signal(signal.SIGTERM,sigstop);signal.signal(signal.SIGINT,sigstop)
if STOP.exists():raise SystemExit('STOP exists: no automatic restart permitted')
if shutil.disk_usage(ROOT).free<sum(int(x['bytes']) for x in PLAN)+1024**3:raise SystemExit('Insufficient conservative disk reserve')
write(ROOT/'process.json',{'pid':os.getpid(),'started':utc(),'proxy':'http://127.0.0.1:7897','max_concurrent':2,'total_bytes':6729410219,'stop_file':str(STOP),'auto_restart':False})
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
 name=item['filename'];size=int(item['bytes']);final=ROOT/name;part=ROOT/(name+'.partial')
 state(name,state='checking_existing',expected_bytes=size,committed_bytes=0,active_chunk_bytes=0)
 try:
  if final.exists():
   if final.stat().st_size==size and md5(final)==item['md5']:
    state(name,state='complete_verified',committed_bytes=size,md5=item['md5']);return
   raise RuntimeError('Existing final file fails size/MD5; preserved without overwrite')
  offset=part.stat().st_size if part.exists() else 0
  if offset>size:raise RuntimeError('Existing partial exceeds expected length')
  # All committed bytes are exact validated contiguous ranges; no test fragments imported.
  while offset<size:
   if stopped():raise InterruptedError('Stop flag or sibling failure')
   end=min(offset+64*1024*1024,size)-1;want=end-offset+1;ok=False
   for attempt in range(1,4):
    if stopped():raise InterruptedError('Stop requested')
    stem=ROOT/'chunks'/f'{name}.{offset}.{attempt}'
    body=Path(str(stem)+'.partial');headers=Path(str(stem)+'.headers');stderr=Path(str(stem)+'.stderr')
    if body.exists():raise RuntimeError('Existing chunk attempt preserved; manual inspection required')
    cmd=['curl','--disable','--proxy','http://127.0.0.1:7897','--noproxy','','--silent','--show-error','--fail','--connect-timeout','20','--max-time','180','--max-filesize',str(want),'--range',f'{offset}-{end}','--dump-header',str(headers),'--output',str(body),item['url']]
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
     with part.open('ab') as dst,body.open('rb') as src:shutil.copyfileobj(src,dst,8*1024*1024);dst.flush();os.fsync(dst.fileno())
     body.unlink();offset+=want;state(name,committed_bytes=offset,active_chunk_bytes=0);ok=True;break
    if code in [401,403,407,429]:raise RuntimeError(f'HTTP {code}; no bypass or retry')
    if p.returncode not in [5,6,7,18,28,35,52,55,56] and code not in [500,502,503,504]:raise RuntimeError('Nontransient response/range/length failure')
    if attempt<3:
     for _ in range(10):
      if stopped():raise InterruptedError('Stop during bounded retry delay')
      time.sleep(.5)
   if not ok:raise RuntimeError('Transient error retry limit reached (3 total attempts)')
  state(name,state='verifying_MD5',committed_bytes=offset,active_chunk_bytes=0)
  digest=md5(part)
  if digest!=item['md5']:raise RuntimeError('Full-file MD5 mismatch; partial preserved')
  if stopped():raise InterruptedError('Stop before final rename')
  os.replace(part,final);state(name,state='complete_verified',md5=digest,committed_bytes=size)
  event({'file':name,'complete_verified':True,'bytes':size,'md5':digest})
 except InterruptedError as e:state(name,state='stopped',error=str(e))
 except Exception as e:
  state(name,state='failed',error=repr(e));event({'file':name,'error':repr(e)});halt.set()
def snapshot(done=False):
 with mutex:d={'utc':utc(),'pid':os.getpid(),'finished':done,'stop_flag':STOP.exists(),'halt':halt.is_set(),'files':json.loads(json.dumps(states)),'child_pids':dict(children)}
 d['verified_files']=sum(v.get('state')=='complete_verified' for v in d['files'].values());d['expected_files']=4;d['expected_bytes']=6729410219
 d['committed_bytes']=sum(v.get('committed_bytes',0) for v in d['files'].values());d['inflight_bytes']=sum(v.get('active_chunk_bytes',0) for v in d['files'].values())
 write(ROOT/'progress.json',d)
(ROOT/'chunks').mkdir(exist_ok=True)
event({'action':'start','pid':os.getpid(),'scope':'e;4 files;R1/R2/I1/I2'})
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 fs=[pool.submit(worker,x) for x in PLAN]
 while not all(f.done() for f in fs):snapshot();time.sleep(5)
 for f in fs:f.result()
snapshot(True);event({'action':'finished','user_stop':STOP.exists()})
