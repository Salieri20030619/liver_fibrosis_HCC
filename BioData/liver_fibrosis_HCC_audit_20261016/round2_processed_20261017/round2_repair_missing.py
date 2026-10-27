import pathlib,json,subprocess,time,datetime,os,re,hashlib,concurrent.futures,shutil
R=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/round2_processed_20261017');W=R/'GSE244832_resume_v2';D=R/'missing_range_repair';D.mkdir(exist_ok=True)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,v):p.write_text(json.dumps(v,indent=2))
state={'pid':os.getpid(),'started_utc':now(),'status':'repairing_missing_ranges','parts':[104,105,119,137,140],'proxy':'disabled','concurrency':2};save(D/'state.json',state)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244832/suppl/GSE244832_hLIVER_processed_files.tar.gz';total=727189398
def get(job):
 idx,a,b=job;p=D/f'{idx}_{a}_{b}.part';log=D/f'{idx}_{a}_{b}.json'
 if p.exists() and log.exists() and json.loads(log.read_text()).get('verified') and p.stat().st_size==b-a+1:return True
 for attempt in range(3):
  tag=f'{idx}_{a}_{b}_attempt{attempt+1}';h=D/(tag+'.headers');err=D/(tag+'.stderr');tmp=D/(tag+'.tmp');start=now()
  with err.open('wb') as e:cp=subprocess.run(['curl','--noproxy','*','--proxy','','-sS','-L','--fail','--connect-timeout','20','--max-time','180','--range',f'{a}-{b}','-D',str(h),'-o',str(tmp),url],stderr=e)
  hs=h.read_text(errors='replace') if h.exists() else '';matches=re.findall(r'Content-Range:\s*bytes\s+(\d+)-(\d+)/(\d+)',hs,re.I);ok=cp.returncode==0 and bool(re.search(r'HTTP/\S+\s+206',hs)) and bool(matches) and tuple(map(int,matches[-1]))==(a,b,total) and tmp.stat().st_size==b-a+1
  rec={'url':url,'range':[a,b],'total':total,'returncode':cp.returncode,'verified':bool(ok),'started_utc':start,'finished_utc':now(),'bytes':tmp.stat().st_size if tmp.exists() else 0,'headers':str(h),'stderr':str(err)};save(D/(tag+'.json'),rec)
  if ok:tmp.replace(p);rec['sha256']=hashlib.sha256(p.read_bytes()).hexdigest();save(log,rec);print(json.dumps(rec),flush=True);return True
  if attempt<2:time.sleep([8,32][attempt])
 return False
try:
 jobs=[];originals={}
 for idx in state['parts']:
  old=json.loads((W/f'{idx:05}.json').read_text());originals[idx]=old;a,b=old['range'];assert not old['verified'] and b-a+1==4194304
  for start in range(a,b+1,1048576):jobs.append((idx,start,min(b,start+1048575)))
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:ok=list(ex.map(get,jobs))
 if not all(ok):state.update(status='still_incomplete',failed_subranges=[list(jobs[i]) for i,v in enumerate(ok) if not v],finished_utc=now());save(D/'state.json',state);raise RuntimeError('Some missing subranges still unavailable; all verified data retained')
 for idx,old in originals.items():
  dest=W/f'{idx:05}.part';assert not dest.exists()
  with dest.open('xb') as f:
   for k,a,b in jobs:
    if k==idx:f.write((D/f'{k}_{a}_{b}.part').read_bytes())
  assert dest.stat().st_size==4194304;shutil.copyfile(W/f'{idx:05}.json',D/f'{idx:05}.original_failure.json');save(W/f'{idx:05}.json',{**old,'verified':True,'bytes':4194304,'returncode':0,'repair_subranges':True,'repair_log_directory':str(D),'finished_utc':now(),'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
 name='GSE244832_hLIVER_processed_files.tar.gz';out=R/name;assert not out.exists();h=hashlib.sha256()
 with out.open('xb') as f:
  for p in [R/(name+'.part')]+[W/f'{i:05}.part' for i in range(171)]:
   if p.parent==W:assert json.loads(p.with_suffix('.json').read_text())['verified']
   with p.open('rb') as src:
    for block in iter(lambda:src.read(8388608),b''):f.write(block);h.update(block)
 assert out.stat().st_size==total
 state.update(status='archive_complete_validating',sha256=h.hexdigest(),bytes=total,updated_utc=now());save(D/'state.json',state)
 save(R/'range_controller_state.json',{'pid':os.getpid(),'status':'repair_complete_validating','utc':now()})
 cp=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python',str(R/'round2_validate_repaired.py')]);assert cp.returncode==0
 cp=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python',str(R/'round2_archive_structural_audit.py')]);assert cp.returncode==0
 state.update(status='completed',finished_utc=now());save(D/'state.json',state)
except Exception as e:
 state.update(error=repr(e),finished_utc=now());save(D/'state.json',state);raise
