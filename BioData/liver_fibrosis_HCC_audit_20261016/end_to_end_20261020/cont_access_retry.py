from pathlib import Path
import requests,json,datetime,hashlib
P=Path(__file__).parent;O=P/'continuation_official_access_retry_v1';O.mkdir(exist_ok=False);s=requests.Session();s.trust_env=False;rows=[]
urls=[('RELA_ensembl','https://rest.ensembl.org/homology/symbol/human/RELA?target_species=mus_musculus&type=orthologues&sequence=none&content-type=application/json'),('RELA_ensembl37','https://grch37.rest.ensembl.org/homology/symbol/human/RELA?target_species=mus_musculus&type=orthologues&sequence=none&content-type=application/json')]
for gene,gid in [('RELA','5970'),('CREB5','9586'),('SOX5','6660')]:
 urls.extend([(gene+'_CTD_detail','https://ctdbase.org/detail.go?type=gene&acc='+gid+'&view=chem'),(gene+'_CTD_query','http://ctdbase.org/basicQuery.go?bqCat=gene&bq='+gene)])
for name,url in urls:
 row={'name':name,'requested_url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  with s.get(url,timeout=(15,30),stream=True) as r:
   content=b''
   for chunk in r.iter_content(65536):
    content+=chunk
    if len(content)>2000000:break
   (O/(name+'.response')).write_bytes(content);row.update(status=r.status_code,final_url=r.url,bytes=len(content),sha256=hashlib.sha256(content).hexdigest(),content_type=r.headers.get('Content-Type'),truncated=len(content)>2000000)
   if 'ensembl' in name and r.status_code==200:
    try:
     j=json.loads(content);row['homologies']=[{'target_id':h['target']['id'],'target_species':h['target']['species'],'type':h['type']} for d in j.get('data',[]) for h in d.get('homologies',[])]
    except Exception as e:row['parse_error']=repr(e)
 except Exception as e:row['raw_error']=repr(e)
 rows.append(row);(O/'attempts.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)
