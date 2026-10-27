import os,sys,json,time,subprocess,selectors,datetime,csv,hashlib
from pathlib import Path
root=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/hammond_speed_20261023')
plan=json.loads(sys.argv[1])
out=root/'bounded_transfer_v1'
out.mkdir(exist_ok=True)
if (out/'started.json').exists():
 print(json.dumps({'status':'already_started_do_not_repeat','path':str(out),'result_exists':(out/'summary.json').exists()}),flush=True);sys.exit(0)
proxy=os.environ.get('https_proxy')
if proxy!='http://127.0.0.1:7897':raise SystemExit('Existing expected HTTP loopback proxy not active; no direct fallback')
cap=8*1024*1024
(out/'started.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'plan':plan,'max_requests':4,'each_seconds':60,'total_seconds':240,'total_bytes':4*cap,'proxy':'http://127.0.0.1:7897','TLS_verify':True},indent=2))
(out/'script.py').write_text(Path(__file__).read_text())
results=[]
for i,item in enumerate(plan,1):
 base=out/f'{i}_{item["filename"]}'
 body=Path(str(base)+'.partial');hdr=Path(str(base)+'.headers');err=Path(str(base)+'.stderr')
 cmd=['curl','--disable','--proxy',proxy,'--noproxy','','--silent','--show-error','--fail','--connect-timeout','15','--max-time','55','--max-filesize',str(cap),'--range',f'0-{cap-1}','--dump-header',str(hdr),'--output','-','--write-out','%{stderr}\nMETRICS\t%{http_code}\t%{http_connect}\t%{time_starttransfer}\t%{time_total}\t%{size_download}\t%{remote_ip}\n',item['url']]
 t=time.monotonic();n=0;first=None;stop=None
 with err.open('wb') as ef,body.open('wb') as bf:
  proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=ef);sel=selectors.DefaultSelector();sel.register(proc.stdout,selectors.EVENT_READ)
  while True:
   if time.monotonic()-t>=59:stop='wall_time_cap';proc.kill();break
   if n==cap:
    try:proc.wait(timeout=min(1,max(.01,59-(time.monotonic()-t))))
    except subprocess.TimeoutExpired:stop='byte_cap';proc.kill()
    break
   ready=sel.select(timeout=min(1,max(.01,59-(time.monotonic()-t))))
   if ready:
    b=os.read(proc.stdout.fileno(),min(65536,cap-n))
    if not b:break
    if first is None:first=time.monotonic()-t
    bf.write(b);n+=len(b)
   elif proc.poll() is not None:break
  if proc.poll() is None:
   try:proc.wait(timeout=max(.01,60-(time.monotonic()-t)))
   except subprocess.TimeoutExpired:proc.kill();stop='wall_time_cap'
  proc.wait();sel.close();proc.stdout.close()
 elapsed=time.monotonic()-t;headers=hdr.read_text(errors='replace') if hdr.exists() else '';stderr=err.read_text(errors='replace');metrics=[l for l in stderr.splitlines() if l.startswith('METRICS\t')]
 fields=metrics[-1].split('\t')[1:] if metrics else []
 ranges=[l.split(':',1)[1].strip() for l in headers.splitlines() if l.lower().startswith('content-range:')]
 code=fields[0] if fields else None;connect=fields[1] if len(fields)>1 else None
 if proc.returncode==0 and n==cap and code=='206' and ranges and ranges[-1].startswith(f'bytes 0-{cap-1}/'):state='success'
 elif connect not in [None,'000','200']:state='proxy_CONNECT_failure'
 elif code in ['403','401','429']:state='HTTP_refused_stop_no_retry'
 elif proc.returncode==7:state='proxy_TCP_failure'
 elif code=='000' and connect=='200':state='post_CONNECT_TLS_or_origin_failure'
 elif code=='200':state='range_ignored_bounded_stop'
 else:state='data_transfer_or_timeout_failure'
 d={'test':i,**item,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'state':state,'curl_returncode':proc.returncode,'http_status':code,'proxy_connect_status':connect,'content_range':ranges[-1] if ranges else None,'received_saved_bytes':n,'elapsed_seconds':elapsed,'first_payload_seconds':first,'curl_time_starttransfer':float(fields[2]) if len(fields)>2 else None,'average_MiB_s':n/1048576/elapsed,'average_Mbps':n*8/1000000/elapsed,'stop_reason':stop,'partial_path':str(body),'partial_SHA256':hashlib.sha256(body.read_bytes()).hexdigest(),'full_fastq_MD5_verified':False,'stderr':stderr}
 (out/f'test_{i}.json').write_text(json.dumps(d,indent=2));results.append(d);print(json.dumps(d),flush=True)
 # A denied target is never retried. All four requests are distinct file URLs.
valid=[r for r in results if r['state']=='success']
s={'results':results,'total_payload_bytes':sum(r['received_saved_bytes'] for r in results),'successful_tests':len(valid),'requests':len(results),'proxy_path':'sgs10:127.0.0.1:7897 -> existing SSH RemoteForward -> Windows:127.0.0.1:7900 (Clash)','VPN_coverage':'unverified; ExpressVPN status alone does not establish process traffic coverage','output':str(out)}
if valid:
 speeds=[r['received_saved_bytes']/r['elapsed_seconds'] for r in valid];pooled=sum(r['received_saved_bytes'] for r in valid)/sum(r['elapsed_seconds'] for r in valid)
 s['estimates']={name:{'bytes':size,'pooled_hours':size/pooled/3600,'fastest_test_hours':size/max(speeds)/3600,'slowest_test_hours':size/min(speeds)/3600} for name,size in [('R1_R2',107177393242),('all_files',126573389087)]}
else:s['estimates']=None
(out/'summary.json').write_text(json.dumps(s,indent=2));print('SUMMARY '+json.dumps(s),flush=True)
