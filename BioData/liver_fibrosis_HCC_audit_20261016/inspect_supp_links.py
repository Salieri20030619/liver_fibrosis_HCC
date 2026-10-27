import pathlib,re,html,json,sys
p=pathlib.Path(sys.argv[1]);s=p.read_text()
for m in re.finditer(r'href="([^"]+(?:\.pdf|\.docx|\.xlsx|\.zip))"',s):
 href=html.unescape(m.group(1))
 if '/bin/' in href:print(json.dumps({'href':href,'context':re.sub('<[^>]+>',' ',s[max(0,m.start()-180):m.end()+400])},ensure_ascii=False))
