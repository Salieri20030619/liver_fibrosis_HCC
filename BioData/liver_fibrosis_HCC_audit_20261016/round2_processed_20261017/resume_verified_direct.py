import pathlib,subprocess,concurrent.futures,hashlib,json,re,sys,datetime,time,math
g,n,total,root,prefix=sys.argv[1:6];total=int(total);root=pathlib.Path(root);prefix=pathlib.Path(prefix);work=root/(g+'_resume_v2');work.mkdir(exist_ok=True);u=sys.argv[6] if len(sys.argv)>6 else f'https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}/suppl/{n}';curl='curl.exe' if sys.platform=='win32' else 'curl'
start=prefix.stat().st_size if prefix.exists() else 0
# The pre-existing prefix is never appended to or overwritten.
chunk=4*1024*1024
jobs=[(i,a,min(total-1,a+chunk-1)) for i,a in enumerate(range(start,total,chunk))]
def get(job):
 i,a,b=job;p=work/f'{i:05}.part';log=work/f'{i:05}.json'
 if log.exists():
  old=json.loads(log.read_text())
  if old.get('verified') and p.exists() and p.stat().st_size==b-a+1:return old
 for attempt in range(3):
  tag=str(time.time_ns());h=work/f'{i:05}.{tag}.headers';e=work/f'{i:05}.{tag}.stderr';tmp=work/f'{i:05}.{tag}.tmp';ts=datetime.datetime.now(datetime.timezone.utc).isoformat()
  actual_url=u
  with e.open('wb') as err:r=subprocess.run([curl,'--noproxy','*','--proxy','','-sS','-L','--fail','--connect-timeout','20','--max-time','600','--range',f'{a}-{b}','-D',str(h),'-o',str(tmp),actual_url],stderr=err)
  hs=h.read_text(errors='replace') if h.exists() else '';ranges=re.findall(r'Content-Range:\s*bytes\s+(\d+)-(\d+)/(\d+)',hs,re.I);ok=r.returncode==0 and bool(re.search(r'HTTP/\S+\s+206',hs)) and ranges and tuple(map(int,ranges[-1]))==(a,b,total) and tmp.stat().st_size==b-a+1
  v={'url':actual_url,'range':[a,b],'total':total,'returncode':r.returncode,'bytes':tmp.stat().st_size if tmp.exists() else 0,'started_utc':ts,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'verified':bool(ok),'response_headers':str(h),'stderr':str(e)}
  if log.exists():(work/f'{i:05}.{tag}.previous.json').write_bytes(log.read_bytes())
  log.write_text(json.dumps(v,indent=2))
  if ok:tmp.replace(p);print(json.dumps({'resource':g,'part':i,'bytes':v['bytes'],'verified':True}),flush=True);return v
  time.sleep(2+attempt*3)
 return v
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:results=list(ex.map(get,jobs))
if all(x['verified'] for x in results):
 out=root/(n+'.complete');hh=hashlib.sha256()
 with out.open('xb') as dest:
  if start:
   with prefix.open('rb') as f:
    for b in iter(lambda:f.read(8*1024*1024),b''):dest.write(b);hh.update(b)
  for i,_,_ in jobs:
   with (work/f'{i:05}.part').open('rb') as f:
    for b in iter(lambda:f.read(8*1024*1024),b''):dest.write(b);hh.update(b)
 assert out.stat().st_size==total
 v={'resource':g,'url':u,'file':str(out),'bytes':total,'sha256':hh.hexdigest(),'prefix_source':str(prefix),'prefix_bytes':start,'range_verified':True,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};(work/'complete.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
else:print(json.dumps({'resource':g,'failed_parts':[i for i,x in enumerate(results) if not x['verified']]}),flush=True)
