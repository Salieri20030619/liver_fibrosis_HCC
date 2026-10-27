from pathlib import Path
import requests,json,hashlib,datetime
import pandas as pd
P=Path(__file__).parent;O=P/'S8_CTD_actual_queries';O.mkdir(exist_ok=False)
genes=pd.read_csv(P/'F8_F9_real_full_cohort_v1/SHAP_importance.tsv',sep='\t').head(3).gene.tolist();rows=[]
s=requests.Session();s.trust_env=False
for gene in genes:
 for host in ['https://ctdbase.org','https://www.ctdbase.org']:
  url=host+'/basicQuery.go?bqCat=gene&bq='+gene
  z={'gene':gene,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
  try:
   with s.get(url,timeout=30,stream=True) as r:
    payload=next(r.iter_content(chunk_size=1000000),b'');z.update(status=r.status_code,content_type=r.headers.get('Content-Type'),bytes_saved=len(payload),sha256=hashlib.sha256(payload).hexdigest(),final_url=r.url)
    (O/(gene+'_'+('www' if 'www.' in host else 'main')+'.html')).write_bytes(payload)
  except Exception as e:z['error']=repr(e)
  rows.append(z);print(json.dumps(z),flush=True)
  if z.get('status')==200:break
pd.DataFrame(rows).to_csv(O/'access_results.tsv',sep='\t',index=False)
(O/'validation.json').write_text(json.dumps({'status':'retrieved_requires_parse' if any(z.get('status')==200 for z in rows) else 'blocked','genes_from_training_chosen_model_SHAP':genes,'records':rows,'limit':'Only targeted official CTD lookup; no drug/causality inference from failed requests; original CTD relationship filter not disclosed'},indent=2))
