import pathlib,json,subprocess,datetime,hashlib,re,csv
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');d=r/'JASPAR_missing_motif_lookup';d.mkdir(exist_ok=False);missing=json.loads((r/'JASPAR_annotation_lookup/lookup_audit.json').read_text())['missing_ids'];logs=[];found=[]
targets=[(m,'https://jaspar.elixir.no/api/v1/matrix/'+m+'/?format=json') for m in missing]
targets.append(('JASPAR2018_CORE_vertebrates_non-redundant.pfm','https://jaspar.elixir.no/download/data/2018/CORE/JASPAR2018_CORE_vertebrates_non-redundant.pfm'))
for name,url in targets:
 p=d/name;ts=datetime.datetime.now(datetime.timezone.utc).isoformat();a=subprocess.run(['curl','-sS','-L','--fail','--connect-timeout','20','--max-time','120','--max-filesize','3000000','-D',str(d/(name+'.headers')),'-o',str(p),url],capture_output=True,text=True);b=p.read_bytes() if p.exists() else b'';logs.append({'url':url,'file':str(p),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'started_utc':ts,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'curl_exit':a.returncode,'stderr':a.stderr})
 if a.returncode:continue
 if name in missing:
  v=json.loads(b)
  if v.get('matrix_id')==name:found.append({'motif_id':name,'name':v.get('name'),'source':url,'tax_group':v.get('tax_group'),'collection':v.get('collection')})
 else:
  for line in b.decode().splitlines():
   if line.startswith('>'):
    parts=line[1:].split(None,1)
    if len(parts)==2 and parts[0] in missing:found.append({'motif_id':parts[0],'name':parts[1],'source':url,'tax_group':'vertebrates (catalog)','collection':'CORE (catalog)'})
v={'missing_before':missing,'found':found,'logs':logs,'still_unresolved':sorted(set(missing)-{x['motif_id'] for x in found}),'limit':'Legacy catalog lookup identifies the exact archived motif version; does not prove which release the GEO authors used.'}
with (d/'lookup_result.json').open('x') as f:json.dump(v,f,indent=2)
print(json.dumps(v),flush=True)
