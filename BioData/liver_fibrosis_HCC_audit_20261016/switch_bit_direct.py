import pathlib,subprocess,os,signal,time,json,datetime
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');d=r/'BIT_range_resume';old=3967592
rows=[s.split(None,5) for s in subprocess.check_output(['ps','-eo','pid,ppid,pgid,sid,stat,args'],text=True).splitlines()[1:]]
group=[x for x in rows if int(x[2])==old]
assert any(int(x[0])==old and str(r/'recover_bit_small_ranges.py') in x[5] for x in group)
assert all(str(r) in x[5] and (int(x[0])==old or x[5].startswith('curl ')) for x in group)
os.killpg(old,signal.SIGTERM)
for _ in range(40):
 live=[x for x in subprocess.check_output(['ps','-eo','pgid,stat'],text=True).splitlines()[1:] if x.split()[0]==str(old) and not x.split()[1].startswith('Z')]
 if not live:break
 time.sleep(.25)
assert not live,'Old process group still active; do not launch replacement'
bad=[]
for p in d.glob('*.part'):
 if not p.name[:5].isdigit():continue
 m=p.with_suffix('.json')
 if not m.exists():bad.append(p.name);continue
 v=json.loads(m.read_text())
 if not v.get('verified') or p.stat().st_size!=v['range'][1]-v['range'][0]+1:bad.append(p.name)
assert not bad,repr(bad)
parts={int(p.stem):json.loads(p.read_text()) for p in d.glob('?????.json')}
a=list(parts.values())+[json.loads(p.read_text()) for p in d.glob('*.sub??.json') if int(p.name[:5]) not in parts]
before=10116144+sum(v['bytes'] for v in a if v.get('verified'))
s=(r/'recover_bit_small_ranges.py').read_text();assert s.count("['curl','-sS'")==1
s=s.replace("['curl','-sS'","['curl','-q','--noproxy','*','-sS'")
s=s.replace("v={'url':endpoint,","v={'transport':'direct; curl -q --noproxy *; proxy environment removed','url':endpoint,")
pre="import os\nfor k in list(os.environ):\n if k.lower() in ('http_proxy','https_proxy','all_proxy','ftp_proxy','no_proxy'):os.environ.pop(k,None)\nos.environ['NO_PROXY']='*'\nos.environ['no_proxy']='*'\n"
script=r/'recover_bit_direct.py'
with script.open('x') as f:f.write(pre+s)
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S');log=(r/f'BIT_direct_{stamp}.log').open('x')
env={k:v for k,v in os.environ.items() if k.lower() not in ('http_proxy','https_proxy','all_proxy','ftp_proxy','no_proxy')};env.update(NO_PROXY='*',no_proxy='*')
p=subprocess.Popen(['python3',str(script)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=env)
v={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'old_pid':old,'old_group_stopped':True,'new_pid':p.pid,'before_verified_bytes':before,'log':log.name,'script':str(script),'proxy_scope':'child process only; no system or SSH settings changed','no_proxy_fallback':True}
with (r/f'BIT_direct_switch_{stamp}.json').open('x') as f:json.dump(v,f,indent=2)
print(json.dumps(v))
