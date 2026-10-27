import pathlib,subprocess,json,datetime,hashlib,os,time,concurrent.futures,gzip,csv,math,tarfile
ROOT=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016')
OUT=ROOT/'round2_processed_20261017';OUT.mkdir(exist_ok=True)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p,v):
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8');tmp.replace(p)
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
TARGETS=[('GSE244832_hLIVER_processed_files.tar.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244832/suppl/GSE244832_hLIVER_processed_files.tar.gz','archive'),('GSE246221_rawcounts_allsamples_human.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246221/suppl/GSE246221_rawcounts_allsamples_human.txt.gz','counts'),('GSE246221_rawcounts_allsamples_mouse.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246221/suppl/GSE246221_rawcounts_allsamples_mouse.txt.gz','counts'),('GSE83596_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83596/matrix/GSE83596_series_matrix.txt.gz','series')]
def validate(p,kind):
 v={'validation_started':now()}
 if kind=='archive':
  members=[];extracted=[];small=OUT/'archive_small_files';small.mkdir(exist_ok=True)
  with tarfile.open(p,'r:gz') as tf:
   for m in tf:
    members.append({'name':m.name,'bytes':m.size,'is_file':m.isfile()})
    if m.isfile() and m.size<=50000000 and any(k in m.name.lower() for k in ['meta','annot','barcode','sample','donor','patient','peak','feature']):
     safe=hashlib.sha256(m.name.encode()).hexdigest()[:8]+'_'+pathlib.PurePosixPath(m.name).name;dest=small/safe
     if not dest.exists():
      with tf.extractfile(m) as src,dest.open('xb') as dst:
       for b in iter(lambda:src.read(1048576),b''):dst.write(b)
     extracted.append({'member':m.name,'path':str(dest),'bytes':dest.stat().st_size,'sha256':digest(dest)})
  save(OUT/(p.name+'.members.json'),members);save(OUT/(p.name+'.extracted_small.json'),extracted)
  # Read the gzip stream through EOF, including CRC/trailer validation.
  total=0
  with gzip.open(p,'rb') as g:
   for b in iter(lambda:g.read(8*1024*1024),b''):total+=len(b)
  v.update(member_count=len(members),extracted_small=len(extracted),gzip_full_EOF_pass=True,uncompressed_stream_bytes=total)
 else:
  n=0;ids=set();dups=widthbad=nonfinite=negative=noninteger=0;headers=None;metadata=[];table=kind=='counts';end=False
  with gzip.open(p,'rt',encoding='utf-8',newline='') as g:
   for line in g:
    if kind=='series' and line.startswith('!'):
     if line.startswith('!series_matrix_table_begin'):table=True
     elif line.startswith('!series_matrix_table_end'):end=True;table=False
     elif line.startswith('!Sample_'):metadata.append(line.rstrip())
     continue
    if not table:continue
    row=next(csv.reader([line],delimiter='\t'))
    if headers is None:headers=row;continue
    if not row:continue
    n+=1;dups+=row[0] in ids;ids.add(row[0]);widthbad+=len(row)!=len(headers)
    for x in row[1:]:
     try:z=float(x)
     except ValueError:nonfinite+=1;continue
     nonfinite+=not math.isfinite(z);negative+=z<0;noninteger+=math.isfinite(z) and z!=int(z)
  cols=headers[1:];save(OUT/(p.name+'.columns.json'),cols)
  if metadata:(OUT/(p.name+'.sample_metadata.txt')).write_text('\n'.join(metadata),encoding='utf-8')
  v.update(shape=[n,len(cols)],row_ID_duplicates=dups,column_ID_duplicates=len(cols)-len(set(cols)),row_width_errors=widthbad,nonfinite=nonfinite,negative=negative,noninteger=noninteger,gzip_full_EOF_pass=True,series_end_marker=end if kind=='series' else None)
  if kind=='counts':
   samples=json.loads((OUT/'GSE246221_SOFT_samples.json').read_text());matches={}
   for s in samples:
    for x in s.get('Sample_description',[]):
     if x in cols:matches.setdefault(x,[]).append(s['GSM'])
   v.update(columns_matching_SOFT=len(matches),unmatched_columns=[x for x in cols if x not in matches],ambiguous_SOFT_columns={k:z for k,z in matches.items() if len(z)!=1})
   save(OUT/(p.name+'.SOFT_column_join.json'),matches)
 v['validation_finished']=now();return v
def work(target):
 name,url,kind=target;p=OUT/name;part=OUT/(name+'.part');state={'filename':name,'url':url,'kind':kind,'started_utc':now(),'proxy':'explicitly disabled','attempts':[],'status':'started'};statepath=OUT/(name+'.state.json');save(statepath,state)
 try:
  if not p.exists():
   for attempt in range(1,5):
    before=part.stat().st_size if part.exists() else 0
    log=OUT/(name+f'.attempt{attempt}.stderr');headers=OUT/(name+f'.attempt{attempt}.headers')
    args=['curl','--noproxy','*','--proxy','','--fail','--location','--connect-timeout','25','--max-time','900','--retry','2','--retry-delay','10','--retry-max-time','180','--continue-at','-','--dump-header',str(headers),'--output',str(part),url]
    with log.open('ab') as f:cp=subprocess.run(args,stdout=f,stderr=f)
    state['attempts'].append({'attempt':attempt,'returncode':cp.returncode,'start_bytes':before,'end_bytes':part.stat().st_size if part.exists() else 0,'finished_utc':now(),'stderr':str(log),'headers':str(headers)});save(statepath,state)
    if cp.returncode==0:part.replace(p);break
    if cp.returncode in [22,33]:break
    if attempt<4:time.sleep(30)
  if not p.exists():state['status']='download_incomplete';return state
  state.update(status='download_complete_validating',bytes=p.stat().st_size,sha256=digest(p),download_finished_utc=now());save(statepath,state)
  state['validation']=validate(p,kind);state['status']='completed'
 except Exception as e:state.update(status='error',error=repr(e))
 finally:state['finished_utc']=now();save(statepath,state)
 return state
if __name__=='__main__':
 save(OUT/'job_identity.json',{'pid':os.getpid(),'started_utc':now(),'root':str(OUT),'targets':[x[0] for x in TARGETS],'proxy':'curl --noproxy * --proxy empty; no network configuration changed'})
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:results=list(ex.map(work,TARGETS))
 save(OUT/'all_download_validation_results.json',results);print(json.dumps(results),flush=True)
