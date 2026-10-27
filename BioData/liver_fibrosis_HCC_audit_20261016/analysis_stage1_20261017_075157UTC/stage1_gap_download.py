from pathlib import Path
import urllib.request,urllib.parse,json,hashlib,gzip,re,datetime,traceback
R=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');O=Path(__file__).parent/'gap_download';O.mkdir(exist_ok=True)
op=urllib.request.build_opener(urllib.request.ProxyHandler({}));records=[]
def save(): (O/'audit.json').write_text(json.dumps(records,indent=2))
def download(url,name):
 p=O/name
 candidates=[q for base in [R,R/'round2_processed_20261017',O] for q in base.glob('*'+name) if q.is_file()]
 for q in candidates:
  if name.endswith('.gz'):
   try:
    with gzip.open(q,'rb') as f:
     while f.read(1048576):pass
    records.append({'url':url,'path':str(q),'status':'reused_CRC_pass','bytes':q.stat().st_size});save();return q
   except Exception:pass
 with op.open(url,timeout=90) as res,open(str(p)+'.partial','wb') as out:
  h=hashlib.sha256();n=0
  while True:
   b=res.read(1048576)
   if not b:break
   out.write(b);h.update(b);n+=len(b)
 if name.endswith('.gz'):
  with gzip.open(str(p)+'.partial','rb') as f:
   while f.read(1048576):pass
 Path(str(p)+'.partial').rename(p)
 records.append({'url':url,'path':str(p),'bytes':n,'sha256':h.hexdigest(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'downloaded_CRC_pass' if name.endswith('.gz') else 'downloaded'});save();return p
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166504/suppl/'
try:
 with op.open(base,timeout=45) as r:html=r.read().decode()
 (O/'GSE166504_listing.html').write_text(html)
 links=re.findall(r'href="([^"]+)"',html)
 for suffix in ['cell_metadata.20220204.tsv.gz','cell_raw_counts.20220204.txt.gz']:
  link=next(x for x in links if x.endswith(suffix));download(urllib.parse.urljoin(base,link),Path(link).name)
except Exception:records.append({'url':base,'status':'error','error':traceback.format_exc()});save()
print(json.dumps(records),flush=True)
