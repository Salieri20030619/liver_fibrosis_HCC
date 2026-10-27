from pathlib import Path
import subprocess,json,hashlib,datetime,tarfile,gzip,io,csv,concurrent.futures
import numpy as np
from scipy.io import mmread
R=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/freeze_priority_20261019');R.mkdir(exist_ok=True)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(name,d):(R/name).write_text(json.dumps(d,indent=2,ensure_ascii=False))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def work(acc,expected):
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{acc}_RAW.tar';p=R/(acc+'_RAW.tar');state={'accession':acc,'url':url,'official_filelist_bytes':expected,'started_utc':now(),'scope':'processed counts archive, not FASTQ'};save(acc+'_state.json',state)
 if not p.exists():
  with (R/(acc+'_curl.stderr.txt')).open('ab') as err:
   z=subprocess.run(['curl','--noproxy','*','-L','--fail','--retry','2','--connect-timeout','20','--max-time','3600','-C','-','-o',str(p)+'.partial',url],stderr=err,stdout=subprocess.DEVNULL)
  state['download_returncode']=z.returncode
  if z.returncode:save(acc+'_state.json',state);return
  Path(str(p)+'.partial').rename(p)
 state.update(bytes=p.stat().st_size,sha256=sha(p),downloaded_utc=now(),official_size_match=p.stat().st_size==expected,status='validating');save(acc+'_state.json',state)
 rows=[]
 with tarfile.open(p) as t:
  members=[m for m in t.getmembers() if m.isfile()];state['archive_files']=len(members);save(acc+'_state.json',state)
  if acc=='GSE289173':
   for m in members:
    if not m.name.endswith('matrix.mtx.gz'):continue
    prefix=m.name[:-len('matrix.mtx.gz')]
    with gzip.GzipFile(fileobj=t.extractfile(m)) as f:a=mmread(f).tocoo()
    def lines(s):
     with gzip.GzipFile(fileobj=t.extractfile(prefix+s)) as f:return f.read().decode().splitlines()
    genes=lines('features.tsv.gz');bars=lines('barcodes.tsv.gz');ids=[s.split('\t')[0] for s in genes]
    rows.append(dict(file=m.name,genes=a.shape[0],cells=a.shape[1],nnz=a.nnz,nonfinite=int((~np.isfinite(a.data)).sum()),negative=int((a.data<0).sum()),noninteger=int((a.data!=np.floor(a.data)).sum()),unique_genes=len(set(ids))==len(ids),unique_barcodes=len(set(bars))==len(bars),dimensions_match=len(ids)==a.shape[0] and len(bars)==a.shape[1],gene_order_sha256=hashlib.sha256('\n'.join(ids).encode()).hexdigest(),gsm=m.name.split('_')[0]));save(acc+'_matrix_rows.json',rows)
  else:
   for m in members:
    if not m.name.endswith('.csv.gz'):continue
    with gzip.GzipFile(fileobj=t.extractfile(m)) as f:
     reader=csv.reader(io.TextIOWrapper(f));head=next(reader);g=[];nr=bad=nonfinite=negative=nonint=0
     for line in reader:
      nr+=1;g.append(line[0]);bad+=len(line)!=len(head)
      try:a=np.asarray(line[1:],dtype=float);nonfinite+=int((~np.isfinite(a)).sum());negative+=int((a<0).sum());nonint+=int((a!=np.floor(a)).sum())
      except:bad+=1
    rows.append(dict(file=m.name,genes=nr,cells=len(head)-1,bad_rows=bad,nonfinite=nonfinite,negative=negative,noninteger=nonint,unique_genes=len(set(g))==len(g),unique_barcodes=len(set(head[1:]))==len(head)-1,gene_order_sha256=hashlib.sha256('\n'.join(g).encode()).hexdigest(),gsm=m.name.split('_')[0],kind='raw' if '-raw_' in m.name else 'normalized'));save(acc+'_matrix_rows.json',rows)
 state.update(status='complete',finished_utc=now(),matrix_files_checked=len(rows),matrix_rows_file=acc+'_matrix_rows.json');save(acc+'_state.json',state)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
 list(ex.map(lambda x:work(*x),[('GSE289173',2076282880),('GSE202379',595322880)]))
