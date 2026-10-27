import pathlib,json,subprocess,time,datetime
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');py='/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python';attempted=set()
def active(pid,needle):
 try:return needle in pathlib.Path(f'/proc/{pid}/cmdline').read_bytes().decode().replace('\0',' ')
 except FileNotFoundError:return False
while True:
 state={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'GSE212837_download_running':active(3195904,'GSE212837_RAW.tar'),'ATAC_backoff_running':active(3340859,'run_atac_backoff.py')}
 for g in ['GSE212837','GSE227265','GSE227265_TF']:
  logs=[]
  for p in (r/(g+'_resume_v2')).glob('?????.json'):
   try:
    value=json.loads(p.read_text());value['part_id']=int(p.stem);logs.append(value)
   except json.JSONDecodeError:pass
  state[g]={'verified_parts':sum(bool(x.get('verified')) for x in logs),'verified_bytes':sum(x.get('bytes',0) for x in logs if x.get('verified')),'failed_parts':sorted(x['part_id'] for x in logs if not x.get('verified')),'download_complete':(r/(g+'_resume_v2')/'complete.json').exists()}
 print(json.dumps(state),flush=True)
 for g,s in [('GSE212837','validate_212837.py'),('GSE227265_TF','validate_tf.py')]:
  if state[g]['download_complete'] and g not in attempted:
   attempted.add(g)
   with (r/(g+'_final_validation_execution.log')).open('x') as log:
    ret=subprocess.run([py,str(r/s)],stdout=log,stderr=subprocess.STDOUT)
   print(json.dumps({'resource':g,'validation_returncode':ret.returncode}),flush=True)
 if not state['GSE212837_download_running'] and not state['ATAC_backoff_running']:
  with (r/'final_download_state.json').open('x') as f:json.dump(state,f,indent=2)
  print('ALL EXISTING DOWNLOAD JOBS REACHED TERMINAL STATE',flush=True);break
 time.sleep(55)
