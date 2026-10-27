import pathlib,tarfile,gzip,json,datetime
root=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016'); p=root/'GSE136103_RAW.tar';rows=[]
try:
 with tarfile.open(p,'r:') as tf:
  for m in tf:
   if not m.isfile():continue
   r={'file':m.name,'declared_bytes':m.size,'archive_status':'partial_download'}
   try:
    f=tf.extractfile(m);gz=gzip.GzipFile(fileobj=f)
    if 'matrix.mtx' in m.name:
     header=None;n=0;maxr=maxc=0
     for b in gz:
      if b.startswith(b'%'):continue
      a=b.split()
      if header is None:header=list(map(int,a));continue
      if len(a)!=3:raise ValueError('bad width')
      rr,cc=int(a[0]),int(a[1]);float(a[2]);maxr=max(maxr,rr);maxc=max(maxc,cc);n+=1
     r.update(dimensions=header,parsed_entries=n,max_row=maxr,max_col=maxc,complete_member=n==header[2],index_bounds_ok=maxr<=header[0] and maxc<=header[1])
    else:
     lines=list(gz);r.update(rows=len(lines),complete_member=True)
   except Exception as e:r['error']=repr(e);rows.append(r);break
   rows.append(r)
except Exception as e:rows.append({'archive_error':repr(e)})
out={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'archive_bytes_at_start':p.stat().st_size,'resource_complete':False,'members':rows};(root/'GSE136103_partial_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
