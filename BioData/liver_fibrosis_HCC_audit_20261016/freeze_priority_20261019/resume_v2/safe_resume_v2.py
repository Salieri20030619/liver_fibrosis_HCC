from pathlib import Path
import subprocess,json,os,time,hashlib,datetime,re,fcntl,concurrent.futures,traceback
R=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/freeze_priority_20261019')
V=R/'resume_v2';V.mkdir(exist_ok=True)
SIZES={'GSE289173':2076282880,'GSE202379':595322880}
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write(path,data):
 tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2));tmp.replace(path)
def header(path):
 blocks=re.split(r'\r?\n\r?\n',path.read_text(errors='replace'))
 blocks=[b for b in blocks if re.match(r'HTTP/\S+\s+\d+',b)]
 if not blocks:return {}
 lines=blocks[-1].splitlines();out={'status':int(lines[0].split()[1])}
 for l in lines[1:]:
  if ':' in l:k,v=l.split(':',1);out[k.lower()]=v.strip()
 return out
def request(acc,start,end,tag,validator=None):
 base=V/(acc+'_'+tag);h=Path(str(base)+'.headers');b=Path(str(base)+'.body');e=Path(str(base)+'.stderr')
 if any(x.exists() for x in [h,b,e]):raise RuntimeError('Refuse overwrite attempt '+str(base))
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{acc}_RAW.tar'
 cmd=['curl','--noproxy','*','--location','--fail','--silent','--show-error','--connect-timeout','30','--max-time','1800','--speed-limit','1','--speed-time','300','--range',f'{start}-{end}','--max-filesize',str(end-start+1),'--dump-header',str(h),'--output',str(b),url]
 if validator:cmd[1:1]=['--header','If-Range: '+validator]
 with e.open('wb') as err:p=subprocess.Popen(cmd,stderr=err,stdout=subprocess.DEVNULL)
 write(V/(acc+'_active_request.json'),{'pid':p.pid,'start':start,'end':end,'body':str(b),'started_utc':now(),'command':cmd})
 rc=p.wait();hh=header(h);n=b.stat().st_size if b.exists() else 0
 cr=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',hh.get('content-range',''))
 valid=hh.get('status')==206 and cr is not None and tuple(map(int,cr.groups()))==(start,end,SIZES[acc]) and 0<n<=end-start+1
 record={'started_offset':start,'requested_end':end,'returncode':rc,'headers':hh,'body_bytes':n,'valid_206_range':valid,'finished_utc':now()};write(Path(str(base)+'.json'),record)
 return b,record
def run(acc):
 lock=(V/(acc+'.lock')).open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 expected=SIZES[acc];part=R/(acc+'_RAW.tar.partial');final=R/(acc+'_RAW.tar');statefile=V/(acc+'_state.json')
 state={'accession':acc,'pid':os.getpid(),'expected_bytes':expected,'started_utc':now(),'state':'checking_existing_prefix','direct':True,'internal_curl_retry':False};write(statefile,state)
 def update(**kw):state.update(kw,updated_utc=now());write(statefile,state)
 if final.exists():raise RuntimeError('Unexpected existing final; inspect before reuse')
 n=part.stat().st_size if part.exists() else 0
 if n>expected:raise RuntimeError('Existing partial larger than official file')
 validator=None
 for label,start,end in [('prefix',0,min(n,4096)-1),('tail',max(0,n-4096),n-1)] if n else [('range_support',0,4095)]:
  b,d=request(acc,start,end,'probe_'+label)
  if not d['valid_206_range'] or d['returncode'] or d['body_bytes']!=end-start+1:raise RuntimeError('Range/206 support probe failed: '+str(d))
  if n:
   with part.open('rb') as f:f.seek(start);existing=f.read(end-start+1)
   if existing!=b.read_bytes():raise RuntimeError('Existing prefix/tail mismatch; original partial preserved')
  v=d['headers'].get('etag') or d['headers'].get('last-modified')
  if validator and v!=validator:raise RuntimeError('Resource changed between probes')
  validator=v
 update(state='downloading_verified_ranges',initial_preserved_bytes=n,committed_bytes=n,validator=validator,existing_prefix_tail_verified=bool(n))
 failures=0;attempt=0
 while n<expected:
  attempt+=1;end=min(n+4*1024*1024,expected)-1
  b,d=request(acc,n,end,f'attempt_{attempt:06d}',validator)
  v=d['headers'].get('etag') or d['headers'].get('last-modified')
  if d['valid_206_range'] and validator and v!=validator:raise RuntimeError('Resource validator changed; response not appended')
  if d['valid_206_range']:
   if part.exists() and part.stat().st_size!=n:raise RuntimeError('Unexpected second writer/length change')
   data=b.read_bytes()
   with part.open('ab') as f:f.write(data);f.flush();os.fsync(f.fileno())
   n+=len(data);failures=0
   # Preserve transferred chunks as evidence; only two approved files, bounded extra disk.
   update(committed_bytes=n,attempt=attempt,last_returncode=d['returncode'],state='downloading_verified_ranges')
  else:
   failures+=1;update(state='retryable_response_not_appended',attempt=attempt,committed_bytes=n,last_response=d,consecutive_failures=failures)
   if d['headers'].get('status')==200:raise RuntimeError('HTTP 200 response rejected; never appended')
   if failures>=3:raise RuntimeError('Three consecutive attempts without validated bytes; partial preserved')
   time.sleep(min(60,10*failures))
 update(state='full_length_archive_check',committed_bytes=n)
 import tarfile
 expected_members={x['Name']:int(x['Size']) for x in json.loads((R/(acc+'_official_filelist.json')).read_text()) if x['#Archive/File']=='File'}
 with tarfile.open(part) as t:actual={x.name:x.size for x in t.getmembers() if x.isfile()}
 if actual!=expected_members:raise RuntimeError('Full-size TAR member names/sizes differ; no matrix analysis launched')
 if part.stat().st_size!=expected:raise RuntimeError('Full file size changed')
 part.rename(final);update(state='matrix_integrity_check')
 # Reuse only previously written integrity functions, without launching their download/wait loops.
 ns={};code=(R/'worker.py').read_text().split('\nwith concurrent.futures.ThreadPoolExecutor')[0];exec(compile(code,'existing_integrity_functions','exec'),ns);ns['work'](acc,expected)
 ns2={};code=(R/'completion_watcher.py').read_text().split("\nsave('completion_watcher_state.json'")[0];exec(compile(code,'existing_clinical_functions','exec'),ns2);result=ns2['check'](acc)
 update(state='complete' if result['file_integrity_pass'] else 'completed_with_validation_failure',final_result=acc+'_final_integrity.json')
def guarded(acc):
 try:run(acc)
 except Exception:
  p=V/(acc+'_state.json');d=json.loads(p.read_text()) if p.exists() else {};d.update(state='blocked',error=traceback.format_exc(),updated_utc=now());write(p,d)
write(V/'launcher.json',{'pid':os.getpid(),'started_utc':now(),'resources':SIZES})
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(guarded,SIZES))
