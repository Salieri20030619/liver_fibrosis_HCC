import pathlib,subprocess,time,json,datetime
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');py='/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python'
def active(pid):
 try:return 'GSE212837_RAW.tar' in pathlib.Path(f'/proc/{pid}/cmdline').read_bytes().decode()
 except FileNotFoundError:return False
while active(3195904):time.sleep(20)
def step(name,args,folder=None):
 with (r/(name+'.log')).open('x') as out:
  p=subprocess.Popen(args,stdout=out,stderr=subprocess.STDOUT)
  while p.poll() is None:
   v={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'step':name,'pid':p.pid,'running':True}
   if folder:
    parts=list((r/folder).glob('?????.part'));v.update(verified_parts=len(parts),verified_bytes=sum(x.stat().st_size for x in parts))
   print(json.dumps(v),flush=True)
   for _ in range(10):
    if p.poll() is not None:break
    time.sleep(5)
 print(json.dumps({'step':name,'returncode':p.returncode}),flush=True)
 return p.returncode
if not (r/'GSE212837_RAW.tar.complete').exists():step('GSE212837_bounded_backoff',['python3',str(r/'resume_backoff.py'),'GSE212837','GSE212837_RAW.tar','1913057280',str(r),str(r/'GSE212837_RAW.absent_prefix')],'GSE212837_resume_v2')
if (r/'GSE212837_RAW.tar.complete').exists() and not (r/'GSE212837_full_validation.json').exists():step('GSE212837_after_backoff_validation',[py,str(r/'validate_212837.py')])
step('BIT_resumable_extraction',['python3',str(r/'bit_resumable_extract.py')],'BIT_range_resume')
if (r/'BIT_range_resume/Liver.h5ad.complete').exists():step('BIT_complete_object_metadata',[py,str(r/'inspect_bit_h5.py')])
v={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'own_downloads_terminal':True,'GSE212837_complete':(r/'GSE212837_RAW.tar.complete').exists(),'GSE212837_parsed':(r/'GSE212837_full_validation.json').exists(),'BIT_complete_object':(r/'BIT_range_resume/Liver.h5ad.complete').exists(),'BIT_metadata_parsed':(r/'BIT_h5_inspection/BIT_h5_metadata_audit.json').exists()}
with (r/'terminal_audit_state.json').open('x') as f:json.dump(v,f,indent=2)
print(json.dumps(v),flush=True)
