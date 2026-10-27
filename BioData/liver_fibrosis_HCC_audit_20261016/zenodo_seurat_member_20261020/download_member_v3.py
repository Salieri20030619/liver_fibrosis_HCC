from pathlib import Path
import subprocess,requests,struct,json,datetime,hashlib,zlib,time,os,fcntl,shutil,concurrent.futures,traceback,threading
ROOT=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/zenodo_seurat_member_20261020')
URL='https://zenodo.org/api/records/14586466/files/Zenodo_250212.zip/content'
TOTAL=6647735352;OFFSET=372495227;COMPRESSED=2424133037;UNCOMPRESSED=2442015643;CRC=0x11738af9
MEMBER='Zenodo/processed_seurat_object/MASLD_snRNA_seq_seurat_v4.rds'
ROOT.mkdir(exist_ok=True);PARTS=ROOT/'parts';PARTS.mkdir(exist_ok=True)
lock=(ROOT/'worker.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,v):
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(v,ensure_ascii=False,indent=2));tmp.replace(p)
def sess():
 s=requests.Session();s.trust_env=False;return s
S=sess();events_lock=threading.Lock()
def event(v):
 with events_lock:
  with (ROOT/'events.jsonl').open('a') as f:f.write(json.dumps({'utc':now(),**v})+'\n')
state={'pid':os.getpid(),'started_utc':now(),'status':'preflight','url':URL,'member':MEMBER,'only_member':True,'direct':True}
def update(**kw):state.update(kw,updated_utc=now());save(ROOT/'state.json',state)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
update(status='preflight_curl_zip_descriptor',transport='curl_direct')
try:
 free=shutil.disk_usage(ROOT).free
 assert free>COMPRESSED+UNCOMPRESSED+2*1024**3,('insufficient_disk',free)
 api=ROOT/'official_record_download.json'
 cp=subprocess.run(['curl','--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','15','--max-time','45','--max-filesize','1000000','--output',str(api),'https://zenodo.org/api/records/14586466'],capture_output=True,text=True)
 assert cp.returncode==0,('record_fetch_failed',cp.returncode,cp.stderr)
 record=json.loads(api.read_text())
 file=next(x for x in record['files'] if x['key']=='Zenodo_250212.zip')
 assert file['size']==TOTAL and file['checksum']=='md5:1fd1a014829a8d6a0f2e97fcb07e69fb'
 save(ROOT/'official_record.json',record)
 h=(ROOT/'curl_probe.body').read_bytes()
 headertext=(ROOT/'curl_probe.headers').read_text()
 headers={l.split(':',1)[0].lower():l.split(':',1)[1].strip() for l in headertext.splitlines() if ':' in l}
 assert '206' in headertext.splitlines()[0] and headers.get('content-range')==f'bytes {OFFSET}-{OFFSET+1023}/{TOTAL}'
 assert len(h)==1024 and h[:4]==b'PK\x03\x04'
 vals=struct.unpack_from('<IHHHHHIIIHH',h);namelen,extra=vals[-2:]
 assert h[30:30+namelen].decode()==MEMBER and vals[3]==8 and vals[8]==UNCOMPRESSED,vals
 if vals[2]&8:
  assert vals[6]==0 and vals[7]==0,vals
  ds=OFFSET+30+namelen+extra+COMPRESSED;dp=ROOT/'member_data_descriptor.bin';dh=ROOT/'member_data_descriptor.headers'
  cp=subprocess.run(['curl','--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','15','--max-time','45','--max-filesize','32','--range',f'{ds}-{ds+31}','--dump-header',str(dh),'--output',str(dp),URL],capture_output=True,text=True)
  assert cp.returncode==0 and f'bytes {ds}-{ds+31}/{TOTAL}' in dh.read_text(),('descriptor_range',cp.returncode,cp.stderr)
  dd=dp.read_bytes();assert len(dd)==32;dv=struct.unpack_from('<IIII',dd)
  assert dv==(0x08074b50,CRC,COMPRESSED,UNCOMPRESSED),dv
 else:assert vals[6]==CRC and vals[7]==COMPRESSED,vals

 START=OFFSET+30+namelen+extra;lowerheaders={k.lower():v for k,v in headers.items()};validator=file['checksum']
 assert validator,('missing_version_validator',headers)
 preflight={'utc':now(),'free_bytes':free,'archive_record_revision':record.get('revision'),'archive_md5':file['checksum'],'headers':headers,'payload_start':START,'payload_end':START+COMPRESSED-1,'compressed_bytes':COMPRESSED,'uncompressed_bytes':UNCOMPRESSED,'crc32':f'{CRC:08x}','validator':validator}
 old=ROOT/'preflight.json'
 if old.exists():assert json.loads(old.read_text())['validator']==validator,'version_changed'
 save(old,preflight);(ROOT/'member_local_header.bin').write_bytes(h)
 CHUNK=8*1024**2;segments=[(i,a,min(COMPRESSED,a+CHUNK)) for i,a in enumerate(range(0,COMPRESSED,CHUNK))]
 update(status='downloading',payload_start=START,compressed_expected=COMPRESSED,segment_count=len(segments),free_bytes=free)
 def fetch(seg):
  i,a,b=seg;final=PARTS/f'{i:04d}.part';meta=PARTS/f'{i:04d}.json';partial=PARTS/f'{i:04d}.partial';size=b-a
  if final.exists():
   m=json.loads(meta.read_text());assert final.stat().st_size==size and sha(final)==m['sha256'];return size
  session=sess()
  for attempt in range(1,5):
   n=partial.stat().st_size if partial.exists() else 0;assert n<=size
   if n==size:break
   lo=START+a+n;hi=START+b-1
   try:
    attemptbase=PARTS/f'{i:04d}_attempt{attempt}_{n}'
    hp=Path(str(attemptbase)+'.headers');bp=Path(str(attemptbase)+'.body')
    cp=subprocess.run(['curl','--noproxy','*','--location','--fail','--silent','--show-error','--connect-timeout','20','--max-time','900','--speed-limit','1024','--speed-time','90','--range',f'{lo}-{hi}','--max-filesize',str(hi-lo+1),'--dump-header',str(hp),'--output',str(bp),URL],capture_output=True,text=True)
    ht=hp.read_text(errors='replace') if hp.exists() else ''
    blocks=[b for b in ht.split('\n\n') if b.startswith('HTTP/')];last=blocks[-1] if blocks else ''
    status=int(last.splitlines()[0].split()[1]) if last else 0
    hdr={l.split(':',1)[0].lower():l.split(':',1)[1].strip() for l in last.splitlines()[1:] if ':' in l}
    received=bp.stat().st_size if bp.exists() else 0
    event({'segment':i,'attempt':attempt,'start':lo,'end':hi,'status':status,'headers':hdr,'returncode':cp.returncode,'stderr':cp.stderr,'received_bytes':received})
    if status==200:raise RuntimeError('FULL_RESPONSE_REFUSED_MAX_FILESIZE_GUARD')
    assert status==206 and hdr.get('content-range')==f'bytes {lo}-{hi}/{TOTAL}',('range_or_version_mismatch',hdr)
    assert 0<=received<=hi-lo+1
    if received:
     with partial.open('ab') as f, bp.open('rb') as src:
      shutil.copyfileobj(src,f,1024*1024);f.flush();os.fsync(f.fileno())
     bp.unlink()
    assert cp.returncode==0 and received==hi-lo+1,('incomplete_valid_prefix_retained',cp.returncode,cp.stderr)
    break
   except Exception as e:
    event({'segment':i,'attempt':attempt,'error':repr(e),'preserved_bytes':partial.stat().st_size if partial.exists() else 0})
    if 'FULL_RESPONSE_REFUSED' in str(e) or 'version_mismatch' in str(e):raise
    if attempt==4:raise
    time.sleep(attempt*3)
  assert partial.stat().st_size==size
  digest=sha(partial);partial.replace(final);save(meta,{'index':i,'relative_start':a,'bytes':size,'absolute_start':START+a,'absolute_end':START+b-1,'sha256':digest,'validator':validator,'completed_utc':now()})
  event({'segment':i,'status':'verified_complete','bytes':size,'sha256':digest});return size
 completed=0
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
  futures={ex.submit(fetch,x):x for x in segments}
  for fut in concurrent.futures.as_completed(futures):
   completed+=fut.result();update(status='downloading',verified_compressed_bytes=completed,verified_segments=sum(1 for p in PARTS.glob('*.part')))
 update(status='decompressing',verified_compressed_bytes=completed)
 out=ROOT/'MASLD_snRNA_seq_seurat_v4.rds';temp=ROOT/'MASLD_snRNA_seq_seurat_v4.rds.partial'
 assert not out.exists(),'existing_final_requires_review'
 dec=zlib.decompressobj(-15);crc=0;total=0;digest=hashlib.sha256();compressed_digest=hashlib.sha256()
 with temp.open('wb') as f:
  for i,a,b in segments:
   with (PARTS/f'{i:04d}.part').open('rb') as src:
    for chunk in iter(lambda:src.read(1024*1024),b''):
     compressed_digest.update(chunk);data=dec.decompress(chunk);assert total+len(data)<=UNCOMPRESSED
     f.write(data);crc=zlib.crc32(data,crc);digest.update(data);total+=len(data)
  tail=dec.flush();f.write(tail);crc=zlib.crc32(tail,crc);digest.update(tail);total+=len(tail);f.flush();os.fsync(f.fileno())
 assert dec.eof and not dec.unused_data and total==UNCOMPRESSED and crc==CRC,(dec.eof,total,hex(crc))
 temp.replace(out);integrity={'completed_utc':now(),'path':str(out),'bytes':total,'crc32':f'{crc:08x}','crc32_matches_zip_directory':True,'sha256':digest.hexdigest(),'compressed_member_sha256':compressed_digest.hexdigest(),'whole_zip_downloaded':False,'whole_zip_md5_verified':False}
 save(ROOT/'member_integrity.json',integrity);update(status='download_and_member_integrity_complete',integrity=integrity)
except Exception as e:
 update(status='blocked_or_failed',error=repr(e),traceback=traceback.format_exc());raise
