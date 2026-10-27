from pathlib import Path
import requests,json,hashlib,gzip,datetime
P=Path(__file__).parent;O=P/'continuation_extra_bulk_metadata_v1';O.mkdir(exist_ok=False);s=requests.Session();s.trust_env=False;rows=[]
for acc in ['GSE167523','GSE126848','GSE66676']:
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz';q=O/(acc+'_family.soft.gz');v={'accession':acc,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};fields=[];sample=''
 try:
  with s.get(url,stream=True,timeout=40) as r:
   v['http_status']=r.status_code;r.raise_for_status();size=0
   with open(q,'wb') as f:
    for b in r.iter_content(65536):
     size+=len(b)
     if size>15000000:raise RuntimeError('Metadata transfer exceeds15MB cap; partialfile notcomplete')
     f.write(b)
  v.update(bytes=q.stat().st_size,sha256=hashlib.sha256(q.read_bytes()).hexdigest());total=0
  with gzip.open(q,'rt',errors='replace') as f:
   for line in f:
    total+=len(line)
    if total>200000000:raise RuntimeError('Decompressed metadata audit exceeds200MB cap')
    if line.startswith('^SAMPLE = '):sample=line.strip().split(' = ',1)[1]
    if line.startswith(('!Series_summary','!Series_overall_design','!Sample_characteristics_ch1','!Sample_title','!Series_supplementary_file','!Sample_supplementary_file')):fields.append((sample,line.strip()))
  v['parse_complete']=True;v['samples']=len(set(x[0] for x in fields if x[0]));v['clinical_keyword_lines']=[x for x in fields if any(k in x[1].lower() for k in ['fibros','cirrho','stage','f0','f1','f2','f3','f4'])];(O/(acc+'_metadata_fields.tsv')).write_text('GSM\tfield\n'+'\n'.join(a+'\t'+b.replace('\t',' ') for a,b in fields)+'\n')
 except Exception as e:v['error']=repr(e);v['parse_complete']=False
 rows.append(v);(O/'requests.json').write_text(json.dumps(rows,indent=2));print(json.dumps(v),flush=True)
