import pathlib,subprocess,json,datetime,hashlib,re
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');d=r/'BIT_resume_probe';d.mkdir(exist_ok=False)
rows=[]
for name in ['BIT_Liver','BIT_Liver_api_retry']:
 p=r/name;x=json.loads((p/'extraction_result.json').read_text());h=(p/'transfer.headers').read_text();e=(p/'transfer.stderr').read_text()
 rows.append({'attempt':name,'url':x['url'],'bytes_received':x['compressed_bytes_consumed_until_target_end'],'object_bytes_written':x['object_bytes_written'],'elapsed_seconds':(datetime.datetime.fromisoformat(x['finished_utc'])-datetime.datetime.fromisoformat(x['started_utc'])).total_seconds(),'configured_max_time_seconds':1800,'client_timeout_reached':False,'curl_exit':18,'stderr':e,'http_headers_without_cookies':'\n'.join(s for s in h.splitlines() if not s.lower().startswith('set-cookie:'))})
url=rows[1]['url'];probes=[]
for label,start in [('prefix',0),('resume',10116144)]:
 end=start+65535;b=d/(label+'.bin');h=d/(label+'.headers');t=datetime.datetime.now(datetime.timezone.utc).isoformat()
 p=subprocess.run(['curl','-sS','-L','--fail','--range',f'{start}-{end}','--connect-timeout','30','--max-time','120','--max-filesize','1048576','-D',str(h),'-o',str(b),url],capture_output=True,text=True)
 heads=h.read_text() if h.exists() else '';body=b.read_bytes() if b.exists() else b'';verified=p.returncode==0 and len(body)==65536 and f'bytes {start}-{end}/2839519229' in heads
 v={'label':label,'start':start,'end':end,'started_utc':t,'curl_exit':p.returncode,'bytes':len(body),'range_verified':verified,'stderr':p.stderr,'headers':heads,'sha256':hashlib.sha256(body).hexdigest()}
 if label=='prefix':
  with (r/'BIT_Liver_api_retry/archive_prefix.tar.gz.partial').open('rb') as f:v['preserved_prefix_64KiB_matches']=body==f.read(65536)
 probes.append(v)
out={'prior_attempts':rows,'probes':probes,'interpretation':'Transfer interruption is acquisition blocker, not absence of repository annotations. Full continuation deferred until GSE212837 ends; no partial H5 parsing.'}
with (d/'resume_evidence.json').open('x') as f:json.dump(out,f,indent=2)
print(json.dumps(out),flush=True)
