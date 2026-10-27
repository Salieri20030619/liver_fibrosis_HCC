"""One-shot file integrity audit; no cell integration, modelling or inference."""
import pathlib,json,time,os,datetime,tarfile,collections,csv,gzip,io,hashlib,math
R=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/round2_processed_20261017')
def save(name,v):
 p=R/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(v,ensure_ascii=False,indent=2));q.replace(p)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
state={'pid':os.getpid(),'started_utc':now(),'status':'waiting_for_archive'};save('structural_audit_state.json',state)
try:
 deadline=time.time()+21600
 while time.time()<deadline:
  c=json.loads((R/'range_controller_state.json').read_text())
  if c['status']=='completed':break
  if c['status']=='range_download_incomplete':raise RuntimeError('Archive download incomplete; verified ranges retained')
  if not (pathlib.Path('/proc')/str(c['pid'])).exists():raise RuntimeError('Download controller exited before completion')
  time.sleep(45)
 else:raise TimeoutError('Archive not ready within six-hour execution limit')
 state.update(status='auditing_matrix_structure',updated_utc=now());save('structural_audit_state.json',state)
 p=R/'GSE244832_hLIVER_processed_files.tar.gz';results=[]
 with tarfile.open(p,'r:gz') as tf:
  for m in tf:
   if not m.isfile() or '.mtx' not in m.name:continue
   raw=tf.extractfile(m);buf=io.BufferedReader(raw);is_gzip=buf.peek(2)[:2]==b'\x1f\x8b';stream=gzip.GzipFile(fileobj=buf) if is_gzip else buf
   txt=io.TextIOWrapper(stream,encoding='utf-8');header=txt.readline().strip();line=txt.readline()
   while line.startswith('%'):line=txt.readline()
   nr,nc,nnz=map(int,line.split());n=bad=negative=nonfinite=duplicates=regress=0;prevcol=0;seen=set();ncol=collections.Counter();nrow=collections.Counter()
   for line in txt:
    if not line.strip():continue
    n+=1
    try:
     rr,cc,v=line.split();ri=int(rr);ci=int(cc);value=float(v)
     bad+=not(1<=ri<=nr and 1<=ci<=nc);negative+=value<0;nonfinite+=not math.isfinite(value)
     if ci!=prevcol:
      regress+=ci<prevcol;seen.clear();prevcol=ci
     duplicates+=ri in seen;seen.add(ri);ncol[ci]+=1;nrow[ri]+=1
    except Exception:bad+=1
    if n%10000000==0:state.update(rows_checked=n,declared_entries=nnz,updated_utc=now());save('structural_audit_state.json',state)
   results.append({'member':m.name,'inner_gzip_magic':is_gzip,'extension_mismatch':m.name.endswith('.gz') and not is_gzip,'header':header,'shape':[nr,nc],'declared_entries':nnz,'observed_entries':n,'entry_count_matches':n==nnz,'invalid_entries':bad,'negative':negative,'nonfinite':nonfinite,'column_order_regressions':regress,'duplicate_coordinates_within_column_runs':duplicates,'global_duplicate_check_complete':regress==0,'nonempty_rows':len(nrow),'nonempty_columns':len(ncol)})
 save('GSE244832_full_matrix_structure.json',results)
 inv=json.loads((R/'complete_small_member_field_audit.json').read_text());metadata=[]
 links=list(csv.DictReader((R/'GSE244832_18_donor_RNA_GSM_F_Figure4C.tsv').open(encoding='utf-8-sig'),delimiter='\t'));rna={x['RNA_ID']:x for x in links}
 for item in inv:
  path=pathlib.Path(item['path'])
  if not item.get('columns'):continue
  op=gzip.open if path.open('rb').read(2)==b'\x1f\x8b' else open
  with op(path,'rt',encoding='utf-8',newline='') as f:
   first=f.readline();delimiter='\t' if '\t' in first else ',';head=next(csv.reader([first],delimiter=delimiter));sets=[set() for _ in head];rna_matches=[collections.Counter() for _ in head];ordered=[hashlib.sha256() for _ in head];n=0
   for row in csv.reader(f,delimiter=delimiter):
    n+=1
    for i,x in enumerate(row[:len(head)]):
     sets[i].add(x)
     ordered[i].update((x+'\n').encode())
     if x in rna:rna_matches[i][x]+=1
  rec={'member':item['member'],'rows':n,'columns':head,'column_unique_counts':dict(zip(head,map(len,sets))),'column_ID_digests':{head[i]:{'ordered':ordered[i].hexdigest(),'set':hashlib.sha256(('\n'.join(sorted(s))+'\n').encode()).hexdigest()} for i,s in enumerate(sets)},'exact_RNA_ID_fields':{head[i]:dict(c) for i,c in enumerate(rna_matches) if c},'matrix_dimension_comparisons':[{'shape':x['shape'],'rows_equal_genes':n==x['shape'][0],'rows_equal_cells':n==x['shape'][1]} for x in results]};metadata.append(rec)
 save('GSE244832_metadata_dimension_RNA_ID_audit.json',metadata)
 state.update(status='completed',finished_utc=now(),matrix_members=len(results));save('structural_audit_state.json',state)
except Exception as e:
 state.update(status='error',error=repr(e),finished_utc=now());save('structural_audit_state.json',state);raise
