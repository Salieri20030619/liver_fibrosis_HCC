import pathlib,re,json,html
p=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/cosmx_mapping_sources/article.html');s=p.read_text();links=[]
for h in re.findall(r'href="([^"]+)"',s):
 if any(x in h.lower() for x in ['links.lww','journals.lww','doi.org/10.1097']):links.append(html.unescape(h))
print(json.dumps(sorted(set(links)),indent=2))
for pat in ['Fig. S9','Figure S9','FOV','four patients','four liver']:
 for m in list(re.finditer(re.escape(pat),s,re.I))[:5]:print(re.sub('<[^>]+>',' ',s[max(0,m.start()-200):m.end()+900]))
