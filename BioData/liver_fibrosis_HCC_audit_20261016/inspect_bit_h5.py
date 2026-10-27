import pathlib,json,hashlib,datetime,csv
import h5py,numpy as np
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');p=r/'BIT_range_resume/Liver.h5ad.complete';assert p.exists() and p.stat().st_size==4829857046
o=r/'BIT_h5_inspection';o.mkdir(exist_ok=False)
def decode(x):return x.decode() if isinstance(x,bytes) else str(x)
def values(x):
 if isinstance(x,h5py.Dataset):return [decode(v) for v in x[:]]
 if 'codes' in x and 'categories' in x:
  cats=[decode(v) for v in x['categories'][:]];return [cats[int(v)] if int(v)>=0 else '' for v in x['codes'][:]]
 return None
with h5py.File(p,'r') as f:
 summary={'file':str(p),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'root_keys':list(f),'matrix_full_values_scanned':False,'metadata':{}}
 for axis in ['obs','var']:
  g=f[axis];key=decode(g.attrs.get('_index','_index'));ids=values(g[key]);cols={k:values(g[k]) for k in g if k!=key};cols={k:v for k,v in cols.items() if v is not None and len(v)==len(ids)}
  summary['metadata'][axis]={'rows':len(ids),'index_key':key,'duplicate_ids':len(ids)-len(set(ids)),'index_order_sha256':hashlib.sha256('\n'.join(ids).encode()).hexdigest(),'index_examples':ids[:10],'columns':{k:{'distinct_count':len(set(v)),'examples':sorted(set(v))[:30]} for k,v in cols.items()}}
  if axis=='obs':
   with (o/'BIT_obs_metadata.tsv').open('x',newline='') as out:
    w=csv.writer(out,delimiter='\t');w.writerow(['cell_id']+list(cols));w.writerows([cell]+[v[i] for v in cols.values()] for i,cell in enumerate(ids))
   tf=list(csv.DictReader((r/'GSE227265_TF_cell_ID_order.tsv').open(),delimiter='\t'));summary['TF_id_table_columns']=list(tf[0]);candidate=[k for k in tf[0] if 'cell' in k.lower() or 'barcode' in k.lower()]
   if len(candidate)==1:
    tfids=[x[candidate[0]] for x in tf];summary['exact_TF_cell_ID_comparison']={'tf_count':len(tfids),'intersection':len(set(tfids)&set(ids)),'same_set':set(tfids)==set(ids),'same_order':tfids==ids}
 for name in ['X','layers']:
  if name not in f:continue
  g=f[name]
  if isinstance(g,h5py.Dataset):summary[name]={'shape':list(g.shape),'dtype':str(g.dtype)}
  else:summary[name]={'attrs':{k:np.asarray(v).tolist() if not isinstance(v,(str,bytes)) else decode(v) for k,v in g.attrs.items()},'keys':list(g),'datasets':{k:{'shape':list(g[k].shape),'dtype':str(g[k].dtype)} for k in g if isinstance(g[k],h5py.Dataset)}}
summary['interpretation']='Exact ID/coordinate comparison required before mapping to original GEO matrices. Reprocessed data is not automatically same peak/order version. Cell suffixes do not establish donors.'
with (o/'BIT_h5_metadata_audit.json').open('x') as out:json.dump(summary,out,indent=2)
print(json.dumps(summary),flush=True)
