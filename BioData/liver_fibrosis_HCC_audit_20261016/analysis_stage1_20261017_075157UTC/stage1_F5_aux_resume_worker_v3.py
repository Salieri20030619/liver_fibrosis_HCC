from pathlib import Path
import requests,json,hashlib,datetime,concurrent.futures,traceback,time
P=Path(__file__).parent;O=P.parent/'reference_cisTarget_aux_v1';src=json.loads((O/'source.json').read_text());u=src['url'];size=src['bytes_expected'];name=u.rsplit('/',1)[1];chunk=4*1024*1024;d=O/'annotation_segments'
try:
 def get(start):
  end=min(size-1,start+chunk-1);dest=d/(str(start)+'.part')
  if dest.exists() and dest.stat().st_size==end-start+1:return
  for attempt in range(3):
   try:
    q=requests.get(u,headers={'Range':f'bytes={start}-{end}'},timeout=(30,180));q.raise_for_status();assert q.status_code==206 and q.headers.get('Content-Range')==f'bytes {start}-{end}/{size}';assert len(q.content)==end-start+1;dest.write_bytes(q.content);return
   except Exception as e:
    with open(O/'resume_range_errors.jsonl','a') as f:f.write(json.dumps({'start':start,'attempt':attempt+1,'error':str(e)})+'\n')
    if getattr(getattr(e,'response',None),'status_code',None)==429 or attempt==2:raise
    time.sleep(2**attempt)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(get,range(0,size,chunk)))
 h=hashlib.sha256()
 with open(O/name,'xb') as f:
  for start in range(0,size,chunk):b=(d/(str(start)+'.part')).read_bytes();assert len(b)==min(chunk,size-start);f.write(b);h.update(b)
 assert (O/name).stat().st_size==size;tf=O/'allTFs_hg38.txt';(O/'complete.json').write_text(json.dumps({'annotation_bytes':size,'annotation_sha256':h.hexdigest(),'TF_bytes':tf.stat().st_size,'TF_sha256':hashlib.sha256(tf.read_bytes()).hexdigest(),'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'official_annotation_checksum':'not provided; ContentRange and size verified, localSHA256 recorded; parse check follows'},indent=2))
except Exception:(O/'error.txt').write_text(traceback.format_exc());raise
