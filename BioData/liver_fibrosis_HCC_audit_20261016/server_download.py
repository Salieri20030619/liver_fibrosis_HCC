import pathlib,subprocess,json,hashlib,datetime,concurrent.futures
ROOT=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016')
jobs=[('GSE212837','GSE212837_all_data_labeled_processed.h5ad.gz'),('GSE136103','GSE136103_RAW.tar')]
def run(item):
 g,n=item;p=ROOT/n;u=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}/suppl/{n}'
 log={'url':u,'filename':n,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 with open(str(p)+'.curl.log','wb') as e:
  r=subprocess.run(['curl','-L','--fail','--retry','2','--connect-timeout','20','--max-time','2700','-D',str(p)+'.headers','-o',str(p),u],stderr=e)
 log.update(returncode=r.returncode,bytes=p.stat().st_size if p.exists() else 0,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 if r.returncode==0:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
  log['sha256']=h.hexdigest()
 (ROOT/(n+'.log.json')).write_text(json.dumps(log,indent=2))
 print(json.dumps(log),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:list(ex.map(run,jobs))
