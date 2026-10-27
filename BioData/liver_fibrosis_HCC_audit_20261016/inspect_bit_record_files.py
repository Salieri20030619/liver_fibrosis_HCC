import pathlib,subprocess,json,hashlib,datetime
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');d=r/'BIT_record_file_audit';d.mkdir(exist_ok=False);url='https://zenodo.org/api/records/14231098';p=d/'record.json';ts=datetime.datetime.now(datetime.timezone.utc).isoformat();a=subprocess.run(['curl','-sS','-L','--fail','--connect-timeout','20','--max-time','90','--max-filesize','3000000','-D',str(d/'record.headers'),'-o',str(p),url],capture_output=True,text=True);b=p.read_bytes() if p.exists() else b'';v={'url':url,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'started_utc':ts,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'returncode':a.returncode,'stderr':a.stderr}
if a.returncode==0:
 x=json.loads(b);v.update(title=x.get('metadata',{}).get('title'),description=x.get('metadata',{}).get('description'),files=x.get('files'),related_identifiers=x.get('metadata',{}).get('related_identifiers'))
with (d/'record_audit.json').open('x') as f:json.dump(v,f,indent=2)
print(json.dumps(v),flush=True)
