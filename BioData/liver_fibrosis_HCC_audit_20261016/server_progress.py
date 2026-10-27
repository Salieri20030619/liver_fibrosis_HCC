import pathlib,json,datetime
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');out={}
for g in ['GSE136103','GSE189175','GSE227265','GSE227265_TF','GSE212837']:
 p=r/(g+'_resume_v2');logs=[]
 for q in p.glob('?????.json'):
  try:logs.append(json.loads(q.read_text()))
  except:pass
 out[g]={'verified_parts':sum(bool(x.get('verified')) for x in logs),'verified_bytes':sum(x['bytes'] for x in logs if x.get('verified')),'failed_or_pending_logs':sum(not x.get('verified',False) for x in logs),'complete':(p/'complete.json').exists()}
out['spatial_files']={p.name:p.stat().st_size for p in r.glob('GSM*') if p.is_file() and not any(k in p.name for k in ['headers','stderr','log.json'])};out['time']=datetime.datetime.now(datetime.timezone.utc).isoformat();print(json.dumps(out))
