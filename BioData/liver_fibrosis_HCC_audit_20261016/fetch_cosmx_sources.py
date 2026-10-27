import pathlib,subprocess,json,datetime,hashlib
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');p=r/'cosmx_mapping_sources';p.mkdir(exist_ok=False)
jobs=[('article.html','https://pmc.ncbi.nlm.nih.gov/articles/PMC12353172/'),('article.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12353172/fullTextXML'),('oa.xml','https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC12353172')]
for n,u in jobs:
 q=p/n;h=p/(n+'.headers');start=datetime.datetime.now(datetime.timezone.utc).isoformat();a=subprocess.run(['curl','-sS','-L','--fail','--connect-timeout','20','--max-time','60','-D',str(h),'-o',str(q),u],capture_output=True)
 v={'url':u,'file':str(q),'started_utc':start,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'returncode':a.returncode,'stderr':a.stderr.decode(errors='replace'),'bytes':q.stat().st_size if q.exists() else 0,'sha256':hashlib.sha256(q.read_bytes()).hexdigest() if q.exists() else ''}
 (p/(n+'.log.json')).write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
