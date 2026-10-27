from pathlib import Path
import requests,json,hashlib,datetime,concurrent.futures,time,traceback,threading,shutil
P=Path(__file__).parent;O=P.parent/'reference_cisTarget_hg38_mc9nr_v3';O.mkdir(exist_ok=False);OLD=P.parent/'reference_cisTarget_hg38_mc9nr_v2';lock=threading.Lock()
base='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/';names=['hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather','hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather'];sizes=[1248772066,1271446522]
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
jobs=[];manifest=[];chunk=8*1024*1024
for name,size in zip(names,sizes):
 shutil.copyfile(OLD/(name+'.sha1sum.txt'),O/(name+'.sha1sum.txt'));d=O/(name+'.segments');d.mkdir();row={'filename':name,'url':base+name,'expected_bytes':size,'chunk_bytes':chunk,'chunks_completed':0,'bytes_received':0,'started_utc':now(),'status':'downloading','official_checksum_type':'SHA1'};manifest.append(row)
 for start in range(0,size,chunk):jobs.append((row,d,start,min(size-1,start+chunk-1)))
def save(): (O/'progress.json').write_text(json.dumps({'utc':now(),'max_connections':16,'files':manifest},indent=2))
save()
def get(job):
 row,d,start,end=job;path=d/(str(start)+'.part')
 for attempt in range(3):
  try:
   r=requests.get(row['url'],headers={'Range':f'bytes={start}-{end}'},timeout=(30,180));
   if r.status_code==429:raise RuntimeError('HTTP429 rate limit; stop rather than bypass')
   r.raise_for_status();assert r.status_code==206 and r.headers.get('Content-Range')==f"bytes {start}-{end}/{row['expected_bytes']}",(r.status_code,r.headers.get('Content-Range'));assert len(r.content)==end-start+1;path.write_bytes(r.content)
   with lock:row['chunks_completed']+=1;row['bytes_received']+=len(r.content);save()
   return
  except Exception as e:
   with lock:
    with open(O/'range_errors.jsonl','a') as f:f.write(json.dumps({'file':row['filename'],'start':start,'end':end,'attempt':attempt,'utc':now(),'error':repr(e)})+'\n')
   if '429' in str(e) or attempt==2:raise
   time.sleep(2*(attempt+1))
try:
 with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:list(pool.map(get,jobs))
 for row in manifest:
  name=row['filename'];h=hashlib.sha256();h1=hashlib.sha1();dest=O/name
  with open(dest,'xb') as f:
   for start in range(0,row['expected_bytes'],chunk):
    data=(O/(name+'.segments')/(str(start)+'.part')).read_bytes();f.write(data);h.update(data);h1.update(data)
  official=(O/(name+'.sha1sum.txt')).read_text().split()[0];assert h1.hexdigest()==official;assert dest.stat().st_size==row['expected_bytes'];row.update(status='complete_official_SHA1_verified',sha1=h1.hexdigest(),sha256=h.hexdigest(),finished_utc=now());save()
 (O/'rankings_complete.json').write_text(json.dumps({'files':manifest,'utc':now()},indent=2))
except Exception:(O/'worker_error.txt').write_text(traceback.format_exc());raise
