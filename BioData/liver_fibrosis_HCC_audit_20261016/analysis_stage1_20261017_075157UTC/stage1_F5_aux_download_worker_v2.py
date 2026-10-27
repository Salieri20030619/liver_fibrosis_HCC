from pathlib import Path
import requests,json,hashlib,datetime,concurrent.futures,shutil,traceback
P=Path(__file__).parent;O=P.parent/'reference_cisTarget_aux_v1';O.mkdir(exist_ok=True)
try:
 tf=P/'F5_regulatory_resource_audit_v1/allTFs_hg38.txt';assert tf.exists();shutil.copyfile(tf,O/'allTFs_hg38.txt')
 name='motifs-v9-nr.hgnc-m0.001-o0.0.tbl';u='https://resources.aertslab.org/cistarget/motif2tf/'+name;r=requests.head(u,timeout=45);r.raise_for_status();size=int(r.headers['Content-Length']);assert 0<size<300000000;chunk=4*1024*1024;d=O/'annotation_segments';d.mkdir();(O/'source.json').write_text(json.dumps({'url':u,'bytes_expected':size,'headers':dict(r.headers),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
 def get(start):
  end=min(size-1,start+chunk-1);q=requests.get(u,headers={'Range':f'bytes={start}-{end}'},timeout=(30,180));q.raise_for_status();assert q.status_code==206 and q.headers.get('Content-Range')==f'bytes {start}-{end}/{size}';assert len(q.content)==end-start+1;(d/(str(start)+'.part')).write_bytes(q.content)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(get,range(0,size,chunk)))
 h=hashlib.sha256()
 with open(O/name,'xb') as f:
  for start in range(0,size,chunk):b=(d/(str(start)+'.part')).read_bytes();f.write(b);h.update(b)
 assert (O/name).stat().st_size==size;(O/'complete.json').write_text(json.dumps({'annotation_bytes':size,'annotation_sha256':h.hexdigest(),'TF_bytes':tf.stat().st_size,'TF_sha256':hashlib.sha256(tf.read_bytes()).hexdigest(),'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'official_annotation_checksum':'not provided; ContentRange and size verified, localSHA256 recorded; parse check follows'},indent=2));print('complete',flush=True)
except Exception:(O/'error.txt').write_text(traceback.format_exc());raise
