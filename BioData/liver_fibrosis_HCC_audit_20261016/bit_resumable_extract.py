import pathlib,subprocess,concurrent.futures,time,datetime,json,hashlib,tarfile,re,shutil
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');d=r/'BIT_range_resume';d.mkdir(exist_ok=True)
prefix=r/'BIT_Liver_api_retry/archive_prefix.tar.gz.partial';total=2839519229;start=prefix.stat().st_size;chunk=4*1024*1024;expected=4829857046;url='https://zenodo.org/api/records/14231098/files/Input_Data.tar.gz/content'
assert start==10116144
assert json.loads((r/'BIT_resume_probe/resume_evidence.json').read_text())['probes'][0]['preserved_prefix_64KiB_matches']
def fetch(i):
 a=start+i*chunk;b=min(total-1,a+chunk-1);log=d/f'{i:05}.json';part=d/f'{i:05}.part'
 if log.exists() and part.exists():
  x=json.loads(log.read_text())
  if x['verified'] and part.stat().st_size==b-a+1:return part
 for k in range(3):
  if k:time.sleep(60*2**(k-1))
  tag=str(time.time_ns());h=d/f'{i:05}.{tag}.headers';tmp=d/f'{i:05}.{tag}.partial';ts=datetime.datetime.now(datetime.timezone.utc).isoformat()
  p=subprocess.run(['curl','-sS','-L','--fail','--connect-timeout','30','--max-time','600','--max-filesize',str(chunk),'--range',f'{a}-{b}','-D',str(h),'-o',str(tmp),url],capture_output=True,text=True)
  hs=h.read_text() if h.exists() else '';n=tmp.stat().st_size if tmp.exists() else 0;ok=p.returncode==0 and n==b-a+1 and f'bytes {a}-{b}/{total}' in hs
  v={'url':url,'range':[a,b],'bytes':n,'verified':ok,'curl_exit':p.returncode,'stderr':p.stderr,'headers':str(h),'started_utc':ts,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
  with (d/f'{i:05}.{tag}.json').open('x') as f:json.dump(v,f,indent=2)
  if ok:
   tmp.rename(part)
   with log.open('x') as f:json.dump(v,f,indent=2)
   print(json.dumps({'part':i,'bytes':n,'verified':True}),flush=True);return part
 raise RuntimeError(f'Bounded retries failed for range {a}-{b}; preserved prior verified parts')
pool=concurrent.futures.ThreadPoolExecutor(max_workers=2);pending={};count=(total-start+chunk-1)//chunk
class Reader:
 def __init__(self):self.current=prefix.open('rb');self.index=0;self.bytes=0
 def read(self,n=-1):
  assert n>=0
  data=bytearray()
  while len(data)<n:
   b=self.current.read(n-len(data))
   if b:data.extend(b);self.bytes+=len(b);continue
   self.current.close()
   if self.index>=count:break
   for j in range(self.index,min(count,self.index+2)):
    if j not in pending:pending[j]=pool.submit(fetch,j)
   path=pending.pop(self.index).result();self.index+=1;self.current=path.open('rb')
  return bytes(data)
reader=Reader();stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%H%M%S');obj=d/f'Liver.{stamp}.h5ad.partial';written=0;sha=hashlib.sha256();success=False;error='';beg=datetime.datetime.now(datetime.timezone.utc).isoformat()
try:
 with tarfile.open(fileobj=reader,mode='r|gz') as t:
  for m in t:
   if m.name.lstrip('./')!='Input_Data/scATAC_seq/Liver.h5ad':continue
   assert m.size==expected
   with t.extractfile(m) as src,obj.open('xb') as dest:
    while True:
     b=src.read(4*1024*1024)
     if not b:break
     dest.write(b);sha.update(b);written+=len(b)
   assert written==expected
   with obj.open('rb') as f:assert f.read(8)==b'\x89HDF\r\n\x1a\n'
   obj.rename(d/'Liver.h5ad.complete');success=True;break
except Exception as e:error=repr(e)
finally:
 reader.current.close();pool.shutdown(wait=True,cancel_futures=True)
v={'url':url,'started_utc':beg,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'object_bytes':written,'expected_object_bytes':expected,'object_sha256':sha.hexdigest() if success else None,'complete_selected_object':success,'compressed_bytes_consumed':reader.bytes,'archive_total_bytes':total,'archive_fully_verified':False,'prefix_reused_bytes':start,'error':error,'scope':'Only complete selected liver H5AD can be parsed; full archive not claimed verified. No unrelated members extracted.'}
with (d/f'result_{stamp}.json').open('x') as f:json.dump(v,f,indent=2)
print(json.dumps(v),flush=True)
