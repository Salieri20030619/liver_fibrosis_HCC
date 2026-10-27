from pathlib import Path
import requests,json,hashlib,datetime,pandas as pd
P=Path(__file__).parent;O=P/'continuation_extra_bulk_inputs_v1';O.mkdir(exist_ok=False);s=requests.Session();s.trust_env=False;rows=[]
files=[('GSE167523','GSE167523_Sample_phenotype_correspondence.xlsx',5000000),('GSE167523','GSE167523_Raw_gene_counts_matrix.txt.gz',50000000),('GSE126848','GSE126848_Gene_counts_raw.txt.gz',50000000)]
for acc,name,cap in files:
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{name}';q=O/name;v={'accession':acc,'url':url,'name':name,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  with s.get(url,stream=True,timeout=45) as r:
   v['http_status']=r.status_code;r.raise_for_status();size=0
   with open(q,'wb') as f:
    for b in r.iter_content(65536):
     size+=len(b)
     if size>cap:raise RuntimeError('Transfer exceeds declaredcap; partialfile notvalidated')
     f.write(b)
  v.update(bytes=q.stat().st_size,sha256=hashlib.sha256(q.read_bytes()).hexdigest(),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),complete_download=True)
  if name.endswith('xlsx'):
   x=pd.read_excel(q,sheet_name=None);v['sheets']={k:{'shape':list(t.shape),'columns':t.columns.astype(str).tolist(),'head':t.head(3).astype(str).to_dict('records')} for k,t in x.items()}
  else:
   x=pd.read_csv(q,sep='\t',nrows=4);v['header']=x.columns.astype(str).tolist();v['head']=x.iloc[:3,:4].astype(str).to_dict('records')
 except Exception as e:v['error']=repr(e)
 rows.append(v);(O/'download_receipts.json').write_text(json.dumps(rows,indent=2));print(json.dumps(v),flush=True)
