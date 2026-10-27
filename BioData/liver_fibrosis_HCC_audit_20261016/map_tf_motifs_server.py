import pathlib,subprocess,json,csv,datetime,hashlib
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');d=r/'JASPAR_annotation_lookup';d.mkdir(exist_ok=False);ids=list(csv.DictReader((r/'GSE227265_TF_motif_row_IDs.tsv').open(),delimiter='\t'))
url='https://jaspar.elixir.no/api/v1/matrix/?collection=CORE&tax_group=vertebrates&page_size=1000&release=2022&format=json';records=[];logs=[]
for page in range(1,6):
 if not url:break
 dest=d/f'page_{page}.json';ts=datetime.datetime.now(datetime.timezone.utc).isoformat();run=subprocess.run(['curl','-sS','-L','--fail','--connect-timeout','20','--max-time','180','--max-filesize','5000000','-D',str(d/f'page_{page}.headers'),'-o',str(dest),url],capture_output=True,text=True)
 b=dest.read_bytes() if dest.exists() else b'';logs.append({'url':url,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'started_utc':ts,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'curl_exit':run.returncode,'stderr':run.stderr})
 if run.returncode:break
 v=json.loads(b);records.extend(v['results']);url=v.get('next')
lookup={x['matrix_id']:x for x in records};rows=[]
for x in ids:
 m=lookup.get(x['motif_id']);rows.append({**x,'matched_exact_version':bool(m),'JASPAR_name':m.get('name','') if m else '', 'annotation_url':m.get('url','') if m else ''})
with (d/'GSE227265_motif_name_exact_join.tsv').open('x',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
v={'source_motif_rows':len(ids),'lookup_records':len(records),'exact_matches':sum(x['matched_exact_version'] for x in rows),'missing_ids':[x['motif_id'] for x in rows if not x['matched_exact_version']],'RUNX_named_annotations':[x for x in rows if 'RUNX' in x['JASPAR_name'].upper()],'fetch_logs':logs,'limits':['MA accessions are not gene symbols; failed lookup is not motif absence','Lookup release 2022 is not proof of original processing release','Exact motif names do not establish donor mapping or mechanism']}
with (d/'lookup_audit.json').open('x') as f:json.dump(v,f,indent=2)
print(json.dumps({k:x for k,x in v.items() if k!='missing_ids'}))
