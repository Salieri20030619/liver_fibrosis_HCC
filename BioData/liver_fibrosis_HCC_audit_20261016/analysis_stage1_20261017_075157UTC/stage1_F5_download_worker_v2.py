from pathlib import Path
import json,datetime,hashlib,concurrent.futures,requests,traceback
P=Path(__file__).parent;O=P.parent/'reference_cisTarget_hg38_mc9nr_v2';O.mkdir(exist_ok=False)
base='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/'
names=['hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather','hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather']
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def get(url,name,expected=None):
 path=O/name;tmp=O/(name+'.partial');assert not path.exists();v={'url':url,'filename':name,'started_utc':now(),'status':'downloading'};audit=O/(name+'.audit.json');audit.write_text(json.dumps(v,indent=2))
 try:
  r=requests.get(url,stream=True,timeout=(30,180));r.raise_for_status();v['headers']=dict(r.headers);n=0;h=hashlib.sha256();h1=hashlib.sha1()
  with open(tmp,'xb') as f:
   for chunk in r.iter_content(1024*1024):
    f.write(chunk);h.update(chunk);h1.update(chunk);n+=len(chunk)
    if n//(32*1024*1024)>(n-len(chunk))//(32*1024*1024):v.update(bytes_received=n,last_update_utc=now());audit.write_text(json.dumps(v,indent=2))
  v.update(bytes=n,sha256=h.hexdigest(),sha1=h1.hexdigest(),finished_utc=now())
  if expected is not None:assert n==expected,(n,expected)
  if name in names:
   matches=[line.split()[0] for line in (O/(name+'.sha1sum.txt')).read_text().splitlines() if line.split() and line.split()[-1].lstrip('*').split('/')[-1]==name];assert len(set(matches))==1,('official checksum not unique',matches);assert h1.hexdigest()==matches[0];v['matches_official_sha1']=True;v['official_sha256_available']=False
  tmp.rename(path);v['status']='downloaded_official_sha1_verified_local_sha256_recorded' if name in names else 'downloaded_local_sha256';audit.write_text(json.dumps(v,indent=2));return v
 except Exception:v.update(status='failed',error=traceback.format_exc());audit.write_text(json.dumps(v,indent=2));raise
try:
 
 for name in names:get(base+name+'.sha1sum.txt',name+'.sha1sum.txt')
 jobs=[(base+names[0],names[0],1248772066),(base+names[1],names[1],1271446522),('https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl','motifs-v9-nr.hgnc-m0.001-o0.0.tbl',None),('https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt','allTFs_hg38.txt',None)]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(lambda t:get(*t),jobs))
 (O/'complete.json').write_text(json.dumps({'status':'complete','utc':now(),'files':results},indent=2));print('complete',flush=True)
except Exception:(O/'worker_error.txt').write_text(traceback.format_exc());raise
