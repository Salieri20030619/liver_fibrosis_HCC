import pathlib,subprocess,datetime,json,os
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016')
procs=subprocess.check_output(['ps','-eo','pid,args'],text=True)
assert not any(str(r)+'/recover_bit_small_ranges.py' in s for s in procs.splitlines())
assert not any(str(r)+'/start_cosmx_independent.py' in s for s in procs.splitlines())
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S')
log=(r/f'BIT_detached_recovery_{stamp}.log').open('x')
p=subprocess.Popen(['python3',str(r/'recover_bit_small_ranges.py')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
script=r/'finish_existing_slide2.py'
source="""import pathlib,time,os,subprocess,json,datetime
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
"""
with script.open('x') as f:f.write(source)
clog=(r/f'CosMx_Slide2_completion_{stamp}.log').open('x')
q=subprocess.Popen(['python3',str(script)],stdin=subprocess.DEVNULL,stdout=clog,stderr=subprocess.STDOUT,start_new_session=True)
v={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'BIT_pid':p.pid,'BIT_log':log.name,'Slide2_existing_download_pid':3624846,'Slide2_completion_only_pid':q.pid,'Slide2_log':clog.name,'scope':'Resume verified BIT parts; no new CosMx downloader. Server log output detached from SSH session.'}
with (r/f'reconnection_recovery_{stamp}.json').open('x') as f:json.dump(v,f,indent=2)
print(json.dumps(v))
