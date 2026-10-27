from pathlib import Path
import urllib.request,re,json,gzip,hashlib,datetime,time,traceback,csv
P=Path(__file__).parent;O=P/'GPL13667_annotation_v3';O.mkdir(exist_ok=False);name='GPL13667_family.soft.gz';f=O/name;part=O/(name+'.partial');url='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13667/soft/'+name;op=urllib.request.build_opener(urllib.request.ProxyHandler({}));events=[];total=None
def record(**kw):
 events.append({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw});(O/'resume_v2_audit.json').write_text(json.dumps(events,indent=2));print(json.dumps(events[-1]),flush=True)
try:
 import shutil
 shutil.copyfile(P/'GPL13667_annotation_v2/GPL13667_family.soft.gz.partial',part)
 assert not f.exists(),'Complete target already exists: do not duplicate'
 for attempt in range(1,7):
  start=part.stat().st_size if part.exists() else 0
  if total and start==total:break
  try:
   req=urllib.request.Request(url,headers={'Range':'bytes='+str(start)+'-','Accept-Encoding':'identity'})
   with op.open(req,timeout=120) as res:
    cr=res.headers.get('Content-Range','');match=re.fullmatch(r'bytes ([0-9]+)-([0-9]+)/([0-9]+)',cr)
    if start>0:assert res.status==206 and match and int(match.group(1))==start,('Cannot safely resume',res.status,cr)
    total=int(match.group(3)) if match else int(res.headers['Content-Length']);record(attempt=attempt,resume_byte=start,expected_total=total,status=res.status)
    with open(part,'ab') as out:
     while True:
      b=res.read(65536)
      if not b:break
      out.write(b)
   assert part.stat().st_size==total,('truncated response',part.stat().st_size,total)
   break
  except Exception:record(attempt=attempt,error=traceback.format_exc(),retained_bytes=part.stat().st_size);time.sleep(10)
 assert total and part.stat().st_size==total,'Transfer incomplete; retained partial for next recovery'
 with gzip.open(part,'rb') as inp:
  while inp.read(1048576):pass
 h=hashlib.sha256()
 with open(part,'rb') as inp:
  while True:
   b=inp.read(1048576)
   if not b:break
   h.update(b)
 part.rename(f);record(status='complete_CRC_pass',bytes=f.stat().st_size,sha256=h.hexdigest(),url=url)

except Exception:record(status='failed',error=traceback.format_exc());raise
