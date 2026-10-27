import pathlib,subprocess,json,datetime,hashlib,re,html,urllib.parse
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');p=r/'official_supplements';p.mkdir(exist_ok=False)
def fetch(url,name):
 q=p/name;start=datetime.datetime.now(datetime.timezone.utc).isoformat();a=subprocess.run(['curl','-sS','-L','--fail','--connect-timeout','20','--max-time','600','--max-filesize','30000000','-D',str(p/(name+'.headers')),'-o',str(q),url],capture_output=True);b=q.read_bytes() if q.exists() else b'';v={'url':url,'file':str(q),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'started_utc':start,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'returncode':a.returncode,'stderr':a.stderr.decode(errors='replace'),'magic_hex':b[:16].hex(),'is_zip':b.startswith(b'PK'),'is_pdf':b.startswith(b'%PDF')};(p/(name+'.log.json')).write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True);return q
s=(r/'cosmx_mapping_sources/article.html').read_text();links=sorted(set(html.unescape(x) for x in re.findall(r'href="([^"]+)"',s) if '/bin/' in x and x.endswith('.docx')))
for link in links:fetch(urllib.parse.urljoin('https://pmc.ncbi.nlm.nih.gov',link),pathlib.PurePosixPath(link).name)
q=fetch('https://pmc.ncbi.nlm.nih.gov/articles/PMC10750269/','ATAC_article.html');s=q.read_text(errors='replace')
links=sorted(set(html.unescape(x) for x in re.findall(r'href="([^"]+)"',s) if '/bin/' in x and re.search(r'\.(pdf|docx|xlsx|zip|txt|r|R|csv|tsv)$',x)))
(p/'ATAC_actual_supplement_hrefs.json').write_text(json.dumps(links,indent=2))
for link in links:fetch(urllib.parse.urljoin('https://pmc.ncbi.nlm.nih.gov',link),pathlib.PurePosixPath(link).name)
