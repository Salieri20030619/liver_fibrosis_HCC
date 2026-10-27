import pathlib,subprocess,tarfile,hashlib,json,datetime,time,io
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');out=r/'BIT_Liver_api_retry';out.mkdir(exist_ok=False)
url='https://zenodo.org/api/records/14231098/files/Input_Data.tar.gz/content';target='./Input_Data/scATAC_seq/Liver.h5ad';start=datetime.datetime.now(datetime.timezone.utc).isoformat();expected=4829857046
class CountReader:
 def __init__(self,f):self.f=f;self.count=0;self.cache=(out/'archive_prefix.tar.gz.partial').open('xb');self.last=time.monotonic()
 def read(self,n=-1):
  b=self.f.read(n);self.count+=len(b);self.cache.write(b);self.cache.flush()
  if time.monotonic()-self.last>45:print(json.dumps({'compressed_bytes_received':self.count,'object_bytes_written':written}),flush=True);self.last=time.monotonic()
  return b
stderr=(out/'transfer.stderr').open('wb');proc=subprocess.Popen(['curl','-sS','-L','--fail','--connect-timeout','30','--max-time','1800','-D',str(out/'transfer.headers'),url],stdout=subprocess.PIPE,stderr=stderr);reader=CountReader(proc.stdout);members=[];success=False;error='';written=0;hh=hashlib.sha256();last=time.monotonic()
try:
 with tarfile.open(fileobj=reader,mode='r|gz') as t:
  for m in t:
   members.append({'name':m.name,'bytes':m.size,'type':str(m.type)});(out/'members_seen.json').write_text(json.dumps(members,indent=2))
   if m.name.lstrip('./')!=target.lstrip('./'):continue
   assert m.isfile() and m.size==expected,(m.name,m.size)
   with t.extractfile(m) as src,(out/'Liver.h5ad.partial').open('xb') as dest:
    while True:
     b=src.read(4*1024*1024)
     if not b:break
     dest.write(b);hh.update(b);written+=len(b)
     if time.monotonic()-last>45:print(json.dumps({'object_bytes_written':written,'expected_object_bytes':expected,'compressed_bytes_consumed':reader.count}),flush=True);last=time.monotonic()
   assert written==expected
   with (out/'Liver.h5ad.partial').open('rb') as f:assert f.read(8)==b'\x89HDF\r\n\x1a\n'
   (out/'Liver.h5ad.partial').rename(out/'Liver.h5ad');success=True;break
except Exception as e:error=repr(e)
finally:
 reader.cache.close();proc.stdout.close()
 if proc.poll() is None:proc.terminate()
 try:code=proc.wait(timeout=15)
 except subprocess.TimeoutExpired:proc.kill();code=proc.wait()
 stderr.close()
v={'url':url,'archive_declared_compressed_bytes':2839519229,'selected_member':target,'object_expected_bytes':expected,'object_bytes_written':written,'compressed_bytes_consumed_until_target_end':reader.count,'object_sha256':hh.hexdigest() if success else None,'extraction_complete':success,'archive_full_integrity_verified':False,'only_selected_member_extracted':True,'members_encountered':members,'curl_returncode_after_early_stream_close':code,'error':error,'started_utc':start,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Read tar prefix through selected Liver.h5ad only; no unrelated member written. Partial stream termination does not validate entire archive.'}
(out/'extraction_result.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
