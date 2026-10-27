from pathlib import Path
import requests,json,hashlib,datetime
P=Path(__file__).parent;O=P/'continuation_extra_bulk_inputs_v1';q=O/'GSE126848_Gene_counts_raw.txt.gz';assert not q.exists();url='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126848&format=file&file=GSE126848%5FGene%5Fcounts%5Fraw%2Etxt%2Egz';v={'url':url,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};s=requests.Session();s.trust_env=False
try:
 with s.get(url,stream=True,timeout=40) as r:
  v['http_status']=r.status_code;r.raise_for_status();sz=0
  with open(q,'wb') as f:
   for b in r.iter_content(65536):
    sz+=len(b)
    if sz>50000000:raise RuntimeError('Response exceeded50MBcap')
    f.write(b)
 v.update(bytes=q.stat().st_size,sha256=hashlib.sha256(q.read_bytes()).hexdigest(),complete_download=True,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());assert q.read_bytes()[:2]==b'\x1f\x8b','Not gzip'
except Exception as e:v['error']=repr(e)
(O/'GSE126848_alternative_receipt.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
