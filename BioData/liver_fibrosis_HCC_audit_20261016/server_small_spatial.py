import json,pathlib,subprocess,hashlib,datetime,concurrent.futures
ROOT=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');spec=json.loads((ROOT/'additional_file_sizes.json').read_text());jobs=[a for a in spec if a['bytes'] and a['bytes']<30000000]
def run(a):
 p=ROOT/a['filename'];log=ROOT/(a['filename']+'.spatial.log.json')
 if log.exists() and json.loads(log.read_text()).get('complete'):return
 tag=datetime.datetime.now(datetime.timezone.utc).strftime('%H%M%S');part=ROOT/(a['filename']+'.'+tag+'.part');h=ROOT/(a['filename']+'.'+tag+'.headers');err=ROOT/(a['filename']+'.'+tag+'.stderr')
 with err.open('wb') as f:r=subprocess.run(['curl','-sS','-L','--fail','--connect-timeout','20','--max-time','1200','-D',str(h),'-o',str(part),a['url']],stderr=f)
 ok=r.returncode==0 and part.exists() and part.stat().st_size==a['bytes']
 if ok:part.replace(p)
 d={'url':a['url'],'filename':str(p),'expected_bytes':a['bytes'],'received_bytes':p.stat().st_size if ok else part.stat().st_size if part.exists() else 0,'complete':ok,'returncode':r.returncode,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest() if ok else ''};log.write_text(json.dumps(d,indent=2));print(json.dumps(d),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:list(ex.map(run,jobs))
