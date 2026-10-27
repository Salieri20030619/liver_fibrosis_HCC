import pathlib,time,os,subprocess,json,datetime
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016')
pid=3624846
while pathlib.Path('/proc/'+str(pid)).exists():
 cmd=pathlib.Path('/proc/'+str(pid)+'/cmdline').read_bytes()
 if b'resume_backoff.py' not in cmd or b'GSM8855043_RDS' not in cmd:break
 time.sleep(60)
p=r/'GSM8855043_Slide2_SeuratObject.RDS.complete'
v={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'complete':p.exists()}
if p.exists() and not (r/'GSM8855043_RDS_inspection').exists():
 v['inspection_returncode']=subprocess.run(['Rscript',str(r/'inspect_cosmx_rds.R'),str(p),'GSM8855043',str(r/'GSM8855043_RDS_inspection')]).returncode
with (r/'Slide2_post_reconnect_state.json').open('x') as f:json.dump(v,f,indent=2)
print(json.dumps(v),flush=True)
