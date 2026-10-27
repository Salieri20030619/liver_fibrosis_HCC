import pathlib,subprocess,json,time,datetime,hashlib
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');d=r/'CosMx_RDS_transport_probe';d.mkdir(exist_ok=False);name='GSM8855042_Slide1_SeuratObject.RDS';rows=[]
urls={'geo':'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8855042&format=file&file='+name,'ftp':'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8855nnn/GSM8855042/suppl/'+name}
for key,url in urls.items():
 h=d/(key+'.headers');b=d/(key+'.bin');beg=time.monotonic();ts=datetime.datetime.now(datetime.timezone.utc).isoformat()
 p=subprocess.run(['curl','-sS','-L','--fail','--range','0-65535','--max-filesize','1048576','--max-time','120','--connect-timeout','30','-D',str(h),'-o',str(b),url],capture_output=True,text=True)
 data=b.read_bytes() if b.exists() else b'';headers=h.read_text() if h.exists() else '';rows.append({'url':url,'started_utc':ts,'seconds':time.monotonic()-beg,'bytes':len(data),'curl_exit':p.returncode,'stderr':p.stderr,'headers':headers,'sha256':hashlib.sha256(data).hexdigest(),'verified_range':len(data)==65536 and 'bytes 0-65535/431908568' in headers})
with (d/'evidence.json').open('x') as f:json.dump(rows,f,indent=2)
print(json.dumps(rows),flush=True)
