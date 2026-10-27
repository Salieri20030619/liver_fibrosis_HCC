import pathlib,subprocess,json,datetime,hashlib,zipfile,re,html
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');d=r/'CosMx_official_alternatives';d.mkdir(exist_ok=False);fn='NIHMS2075194-supplement-Supplemental_Digital_Content.docx'
jobs=[('europepmc_bin.docx','https://europepmc.org/articles/PMC12353172/bin/'+fn),('europepmc_supplementaryFiles.zip','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12353172/supplementaryFiles'),('pmc_standard_bin.docx','https://pmc.ncbi.nlm.nih.gov/articles/PMC12353172/bin/'+fn),('publisher.html','https://journals.lww.com/10.1097/HEP.0000000000001360')]
for name,url in jobs:
 p=d/name;start=datetime.datetime.now(datetime.timezone.utc).isoformat();a=subprocess.run(['curl','-sS','-L','--fail','--connect-timeout','20','--max-time','240','--max-filesize','30000000','-D',str(d/(name+'.headers')),'-o',str(p),url],capture_output=True,text=True);b=p.read_bytes() if p.exists() else b'';v={'url':url,'file':str(p),'started_utc':start,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'returncode':a.returncode,'stderr':a.stderr,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'magic_hex':b[:16].hex(),'is_zip':zipfile.is_zipfile(p) if p.exists() else False}
 if v['is_zip']:
  with zipfile.ZipFile(p) as z:v['zip_members']=z.namelist();v['zip_crc_pass']=z.testzip() is None;v['is_docx']='word/document.xml' in z.namelist()
 if name=='publisher.html' and b:
  s=b.decode(errors='replace');v['supplement_links']=[{'url':html.unescape(m.group(1)),'context':re.sub('<[^>]+>',' ',s[max(0,m.start()-200):m.end()+200])} for m in re.finditer(r'href=[\"\x27]([^\"\x27]+)[\"\x27]',s) if 'links.lww.com' in m.group(1) or 'supplement' in m.group(1).lower()]
 with (d/(name+'.log.json')).open('x') as f:json.dump(v,f,indent=2)
 print(json.dumps(v),flush=True)
