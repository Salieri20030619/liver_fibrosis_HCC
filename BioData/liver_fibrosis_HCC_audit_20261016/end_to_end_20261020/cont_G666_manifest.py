from pathlib import Path
import requests,json,datetime,hashlib
P=Path(__file__).parent;O=P/'continuation_G666_public_inventory_v1';O.mkdir(exist_ok=False);s=requests.Session();s.trust_env=False;rows=[]
for label,url in [('series','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66676&targ=self&form=text&view=full'),('matrix_listing','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66676/matrix/'),('supp_listing','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66676/suppl/')]:
 v={'label':label,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  r=s.get(url,timeout=40);v['HTTP_status']=r.status_code;r.raise_for_status();assert len(r.content)<2000000;(O/(label+'.txt')).write_bytes(r.content);v.update(bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
 except Exception as e:v['error']=repr(e)
 rows.append(v);(O/'requests.json').write_text(json.dumps(rows,indent=2));print(json.dumps(v),flush=True)
