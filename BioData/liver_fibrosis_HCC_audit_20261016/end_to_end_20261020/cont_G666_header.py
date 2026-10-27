from pathlib import Path
import requests,zlib,json,hashlib,datetime
P=Path(__file__).parent;O=P/'continuation_G666_matrix_header_v1';O.mkdir(exist_ok=False)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66676/matrix/GSE66676_series_matrix.txt.gz'
v={'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'complete_expression_download':False};s=requests.Session();s.trust_env=False;raw=b'';text=b''
try:
 with s.get(url,stream=True,timeout=40) as r:
  v['HTTP_status']=r.status_code;r.raise_for_status();v['content_length']=r.headers.get('Content-Length');d=zlib.decompressobj(16+zlib.MAX_WBITS)
  for b in r.iter_content(16384):
   raw+=b;text+=d.decompress(b)
   if b'!series_matrix_table_begin' in text:break
   if len(raw)>524288:raise RuntimeError('Header-only compressed cap512KiB reached')
 assert b'!series_matrix_table_begin' in text
 header=text.split(b'!series_matrix_table_begin')[0];(O/'header.txt').write_bytes(header)
 fields=[]
 for line in header.decode().splitlines():
  if line.startswith(('!Sample_geo_accession','!Sample_title','!Sample_characteristics_ch1')):
   fields.append(line)
 (O/'sample_fields.tsv').write_text('\n'.join(fields)+'\n')
 v.update(header_complete=True,header_bytes=len(header),header_sha256=hashlib.sha256(header).hexdigest(),sample_fields=fields)
except Exception as e:v.update(header_complete=False,error=repr(e))
(O/'compressed_prefix.bin').write_bytes(raw);v.update(prefix_bytes=len(raw),prefix_sha256=hashlib.sha256(raw).hexdigest(),scope='Only compressed prefix/header inspected; full expression and whole-file checksum NOT verified.')
(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
