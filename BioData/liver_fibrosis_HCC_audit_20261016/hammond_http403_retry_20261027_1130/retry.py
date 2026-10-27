"""Transport only: max three attempts per block, 60s reads, 10/30s backoff; no refusal retries."""
import os,json,time,datetime,fcntl,signal,hashlib,requests,threading
from pathlib import Path
from verified_chunk_io import append_verified
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
lock=(ROOT/'download.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
shared=(BASE/'hammond_e_download_20261023/download.lock').open('a');fcntl.flock(shared,fcntl.LOCK_EX|fcntl.LOCK_NB)
PLAN=json.loads((ROOT/'manifest.json').read_text());assert len(PLAN)==5
STOP=BASE/'hammond_overnight_20261024/STOP';OTHER=BASE/'hammond_e_download_20261023/STOP'
halt=threading.Event();mutex=threading.Lock()
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def stopped():return halt.is_set() or STOP.exists() or OTHER.exists()
def stop(sig,frame):halt.set();STOP.touch()
signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
assert not stopped()
states={x['filename']:{'state':'queued','initial_bytes':x['initial_bytes'],'committed_bytes':x['initial_bytes'],'expected_bytes':int(x['bytes']),'requests':0} for x in PLAN}
def update(name,**kw):
 with mutex:states[name].update(kw)
def event(**kw):
 with (ROOT/'events.jsonl').open('a') as f:f.write(json.dumps({'utc':utc(),**kw})+'\n')
def snapshot(done=False):
 with mutex:d={'utc':utc(),'pid':os.getpid(),'finished':done,'STOP':STOP.exists(),'other_STOP':OTHER.exists(),'halt':halt.is_set(),'files':json.loads(json.dumps(states))}
 d['verified_files']=sum(v['state']=='complete_verified' for v in d['files'].values());d['initial_bytes']=sum(x['initial_bytes'] for x in PLAN)
 d['committed_bytes']=sum(v['committed_bytes'] for v in d['files'].values());d['new_committed_bytes']=d['committed_bytes']-d['initial_bytes']
 t=ROOT/'progress.tmp';t.write_text(json.dumps(d,indent=2));os.replace(t,ROOT/'progress.json')
def ticker():
 while not tickstop.wait(5):snapshot()
tickstop=threading.Event();thread=threading.Thread(target=ticker,daemon=True)
(ROOT/'chunks').mkdir(exist_ok=True)
session=requests.Session();session.trust_env=False;session.proxies={'http':'http://127.0.0.1:17898','https':'http://127.0.0.1:17898'}
session.headers.update({'Accept-Encoding':'identity'})
snapshot();thread.start()
for item in PLAN:
 name=item['filename'];size=int(item['bytes']);final=Path(item['destination'])/name;part=Path(str(final)+'.partial');offset=item['initial_bytes']
 if stopped():update(name,state='stopped');continue
 try:
  assert not final.exists() and (part.stat().st_size if part.exists() else 0)==offset
  while offset<size:
   if stopped():raise InterruptedError('User STOP; no restart')
   for attempt in range(1,4):
    try:
     end=min(size,offset+1048576)-1;want=end-offset+1;stem=ROOT/'chunks'/f'{name}.{offset}.attempt{attempt}'
     body=Path(str(stem)+'.partial');assert not body.exists()
     update(name,state='downloading',range_start=offset,range_end=end,active_bytes=0,requests=states[name]['requests']+1)
     started=time.monotonic()
     with session.get(item['url'],headers={'Range':f'bytes={offset}-{end}'},stream=True,timeout=(20,60),allow_redirects=False,verify=True) as r:
      headers=dict(r.headers);Path(str(stem)+'.headers.json').write_text(json.dumps({'status':r.status_code,'headers':headers},indent=2))
      if r.status_code!=206:
       evidence=b''
       for data in r.iter_content(1024):
        if stopped() or time.monotonic()-started>30:break
        evidence+=data[:8192-len(evidence)]
        if len(evidence)>=8192:break
       ctype=r.headers.get('Content-Type','')
       if any(x in ctype.lower() for x in ['text','html','json','xml']):Path(str(stem)+'.error.txt').write_text(evidence.decode('utf-8',errors='replace'))
       event(file=name,http=r.status_code,range=[offset,end],accepted=False,error_body_bytes=len(evidence),headers=str(stem)+'.headers.json')
       update(name,state='refused' if r.status_code in [401,403,407,429] else 'failed',http=r.status_code,error=f'HTTP {r.status_code}; no retry or redirect',active_bytes=0)
       break
      expected=f'bytes {offset}-{end}/{size}'
      assert r.headers.get('Content-Range')==expected,'Content-Range mismatch'
      assert r.headers.get('Content-Encoding','identity')=='identity','Unexpected encoding'
      assert not r.headers.get('Content-Length') or int(r.headers['Content-Length'])==want
      count=0
      with body.open('xb') as f:
       for data in r.iter_content(65536):
        if stopped():raise InterruptedError('User STOP')
        if time.monotonic()-started>180:raise TimeoutError('Bounded request deadline exceeded')
        count+=len(data);assert count<=want,'Body exceeds requested range';f.write(data);update(name,active_bytes=count)
      assert count==want,'Truncated response'
     proof=append_verified(body,part,offset,want);event(file=name,http=206,range=[offset,end],accepted=True,**proof)
     offset=end+1;update(name,committed_bytes=offset,active_bytes=0);body.unlink()
     break
    except (requests.exceptions.Timeout,requests.exceptions.ConnectionError,requests.exceptions.ChunkedEncodingError,TimeoutError) as e:
     event(file=name,range=[offset,end],attempt=attempt,accepted=False,transport_error=repr(e),retry=attempt<3)
     update(name,active_bytes=0,last_transport_error=repr(e),block_attempt=attempt)
     if attempt==3:raise
     for _ in range(10 if attempt==1 else 30):
      if stopped():raise InterruptedError('User STOP during backoff')
      time.sleep(1)
   if states[name]['state'] in ['refused','failed']:break
  if offset==size:
   update(name,state='verifying_md5');h=hashlib.md5()
   with part.open('rb') as f:
    while True:
     if stopped():raise InterruptedError('User STOP during MD5')
     chunk=f.read(8388608)
     if not chunk:break
     h.update(chunk)
   digest=h.hexdigest();event(file=name,action='full_md5',expected=item['md5'],observed=digest,bytes=part.stat().st_size)
   assert digest==item['md5'] and part.stat().st_size==size,'Full file MD5/length mismatch'
   if stopped():raise InterruptedError('User STOP before rename')
   os.replace(part,final);update(name,state='complete_verified',md5=digest)
 except Exception as e:
  update(name,state='stopped' if isinstance(e,InterruptedError) else 'failed',error=repr(e),active_bytes=0);event(file=name,error=repr(e),accepted=False)
 snapshot()
tickstop.set();thread.join();session.close();snapshot(True);event(action='finished')
shared.close();lock.close()
