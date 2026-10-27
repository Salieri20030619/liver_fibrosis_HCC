from pathlib import Path
import urllib.request,json,hashlib,datetime,traceback,gzip
O=Path(__file__).parent/'GPL13667_annotation_v2';O.mkdir(exist_ok=False)
url='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13667/soft/GPL13667_family.soft.gz'
try:
 op=urllib.request.build_opener(urllib.request.ProxyHandler({}));p=O/'GPL13667_family.soft.gz.partial';h=hashlib.sha256();n=0
 with op.open(url,timeout=120) as r,open(p,'wb') as f:
  while True:
   b=r.read(1048576)
   if not b:break
   f.write(b);h.update(b);n+=len(b);print(json.dumps({'bytes':n}),flush=True)
 text=gzip.open(p,'rt').read();lines=text.splitlines();start=lines.index('!platform_table_begin')+1;end=lines.index('!platform_table_end');rows=lines[start:end];headers=rows[0].split('\t');ids=[r.split('\t')[0] for r in rows[1:]]
 assert len(ids)==49386 and len(set(ids))==49386,(len(ids),len(set(ids)))
 p.rename(O/'GPL13667_family.soft.gz');(O/'probe_annotation.tsv').write_text('\n'.join(rows)+'\n')
 (O/'result.json').write_text(json.dumps({'url':url,'bytes':n,'sha256':h.hexdigest(),'probes':len(ids),'unique_probes':len(set(ids)),'columns':headers,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
except Exception:(O/'error.txt').write_text(traceback.format_exc());raise
