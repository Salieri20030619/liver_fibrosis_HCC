import pathlib,re,html,json
p=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/cosmx_mapping_sources/article.html');s=p.read_text();plain=html.unescape(re.sub('<[^>]+>',' ',s));title=re.search('<title>(.*?)</title>',s,re.S)
links=[]
for x in re.finditer(r'href="([^"]+)"',s):
 u=html.unescape(x.group(1))
 if any(k in u.lower() for k in ['github','zenodo','figshare','mendeley','osf.io','data','geo/query']):links.append({'url':u,'context':html.unescape(re.sub('<[^>]+>',' ',s[max(0,x.start()-100):x.end()+200]))})
contexts=[]
for x in re.finditer('data availability|code availability|accession|github|GSE292268',plain,re.I):contexts.append(plain[max(0,x.start()-200):x.end()+800])
print(json.dumps({'title':title.group(1) if title else '', 'links':links,'contexts':contexts},ensure_ascii=False))
