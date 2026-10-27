import pathlib,subprocess,time,json,datetime
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016')
while not (r/'terminal_audit_state.json').exists():time.sleep(20)
results=[]
for gsm,slide,total in [('GSM8855042','Slide1',431908568),('GSM8855043','Slide2',425519255)]:
 name=f'{gsm}_{slide}_SeuratObject.RDS';url=f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8855nnn/{gsm}/suppl/{name}';g=gsm+'_RDS';prefix=r/'CosMx_RDS_transport_probe/ftp.bin' if slide=='Slide1' else r/'CosMx_Slide2_absent_prefix'
 with (r/(g+'_download_execution.log')).open('x') as log:
  p=subprocess.Popen(['python3',str(r/'resume_backoff.py'),g,name,str(total),str(r),str(prefix),url],stdout=log,stderr=subprocess.STDOUT)
  while p.poll() is None:
   parts=list((r/(g+'_resume_v2')).glob('?????.part'));print(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'resource':g,'running':True,'verified_parts':len(parts),'verified_bytes':sum(x.stat().st_size for x in parts)}),flush=True)
   for _ in range(10):
    if p.poll() is not None:break
    time.sleep(5)
 result={'resource':g,'complete':(r/(name+'.complete')).exists(),'download_returncode':p.returncode}
 if result['complete']:
  with (r/(g+'_R_inspection.log')).open('x') as log:ret=subprocess.run(['Rscript',str(r/'inspect_cosmx_rds.R'),str(r/(name+'.complete')),gsm,str(r/(g+'_inspection'))],stdout=log,stderr=subprocess.STDOUT)
  result['R_returncode']=ret.returncode
 results.append(result);print(json.dumps(result),flush=True)
with (r/'CosMx_RDS_terminal_state.json').open('x') as f:json.dump({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':results},f,indent=2)
