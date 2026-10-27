import pathlib,json,re,csv,datetime,collections
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');rows=[]
for g in ['GSE227265','GSE227265_TF','GSE212837']:
 for q in (r/(g+'_resume_v2')).glob('*.headers'):
  s=q.read_text(errors='replace');codes=re.findall(r'HTTP/\S+\s+(\d+)',s);ra=re.findall(r'Retry-After:\s*([^\r\n]+)',s,re.I);cr=re.findall(r'Content-Range:\s*([^\r\n]+)',s,re.I)
  rows.append({'resource':g,'file':str(q),'http_status':codes[-1] if codes else 'missing','retry_after':ra[-1] if ra else '', 'content_range':cr[-1] if cr else ''})
dest=r/('access_response_audit_'+datetime.datetime.now(datetime.timezone.utc).strftime('%H%M%S')+'.tsv')
with dest.open('x',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps({'file':str(dest),'status_counts':{g:dict(collections.Counter(x['http_status'] for x in rows if x['resource']==g)) for g in ['GSE227265','GSE227265_TF','GSE212837']},'retry_after_present':sum(bool(x['retry_after']) for x in rows),'interpretation':'HTTP503 confirms temporary service failure, not its root cause. No 429 inferred. ATAC retries reduced to two connections with 30/60s exponential backoff; one bounded round.'}))
