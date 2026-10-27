from pathlib import Path
import subprocess,json,tarfile,gzip,csv,collections,datetime
P=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/zenodo_seurat_member_20261020')
O=P/'metadata_audit'
rscript='''md<-readRDS("%s")
z<-gzfile("%s","wt")
write.table(data.frame(cell_id=rownames(md),DonorID=md$DonorID,Disease_group=md$Disease_group),z,sep="\\t",quote=FALSE,row.names=FALSE)
close(z)
write.csv(as.data.frame(table(md$DonorID)),"%s",row.names=FALSE)
'''%(O/'metadata_only.rds',O/'cell_donor_ids.tsv.gz',O/'donor_cell_counts.csv')
r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/r_seurat/bin/Rscript','--vanilla','-'],input=rscript,text=True,capture_output=True);assert r.returncode==0,r.stderr
mapping=json.loads((P/'GEO_donor_map.json').read_text());by={x['donor']:x for x in mapping};cells=collections.defaultdict(list);groups=collections.defaultdict(set)
with gzip.open(O/'cell_donor_ids.tsv.gz','rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):cells[x['DonorID']].append(x['cell_id']);groups[x['DonorID']].add(x['Disease_group'])
rows=[]
with tarfile.open('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/freeze_priority_20261019/GSE289173_RAW.tar','r:') as tar:
 for donor,ids in cells.items():
  m=by.get(donor);prefix=donor+'_'+donor+'_';bad=[x for x in ids if not x.startswith(prefix)];decoded=[x[len(prefix):] for x in ids if x.startswith(prefix)]
  if m:
   name=m['GSM']+'_'+donor+'_barcodes.tsv.gz';raw=gzip.decompress(tar.extractfile(name).read()).decode().splitlines();rawset=set(raw);missing=[x for x in decoded if x not in rawset]
  else:name='';raw=[];missing=decoded
  rows.append({'DonorID':donor,'GSM':m['GSM'] if m else '', 'metadata_cells':len(ids),'raw_cells':len(raw),'raw_barcode_member':name,'exact_prefix':prefix,'unrecognized_prefix_count':len(bad),'decoded_unique_barcodes':len(set(decoded)),'decoded_barcodes_missing_from_raw':len(missing),'metadata_groups':';'.join(sorted(groups[donor])),'GEO_group':m['disease group'] if m else '', 'identity_pass':bool(m) and not bad and not missing and len(set(decoded))==len(ids)})
with (O/'donor_GSM_barcode_identity.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'metadata_donors':len(cells),'GEO_donors':len(by),'donor_sets_equal':set(cells)==set(by),'metadata_cells':sum(len(x) for x in cells.values()),'raw_cells':sum(x['raw_cells'] for x in rows),'all_donor_identity_pass':all(x['identity_pass'] for x in rows),'rule':'Require exact DonorID_DonorID_ prefix, then exact raw barcode membership for the GSM_DonorID named barcode file. No fuzzy matching or aggregate inference.'}
(O/'barcode_identity_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
