from pathlib import Path
import csv,json,gzip,tarfile,re,collections,datetime
P=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/zenodo_seurat_member_20261020');O=P/'barcode_followup';O.mkdir(exist_ok=True)
mapping=json.loads((P/'GEO_donor_map.json').read_text());by={x['donor']:x for x in mapping};raw={};members=[]
with tarfile.open('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/freeze_priority_20261019/GSE289173_RAW.tar','r:') as t:
 for donor,x in by.items():
  name=x['GSM']+'_'+donor+'_barcodes.tsv.gz';m=t.getmember(name);lines=gzip.decompress(t.extractfile(m).read()).decode().splitlines();raw[donor]=set(lines);members.append({'donor':donor,'GSM':x['GSM'],'member':name,'bytes_gzip':m.size,'tar_mtime_utc':datetime.datetime.fromtimestamp(m.mtime,datetime.timezone.utc).isoformat(),'barcode_count':len(lines),'unique_count':len(set(lines)),'suffixes':';'.join(sorted(set(re.sub(r'^[ACGT]+','',s) for s in lines)))})
index=collections.defaultdict(set);core=collections.defaultdict(set)
for donor,barcodes in raw.items():
 for b in barcodes:index[b].add(donor);core[re.sub(r'-\d+$','',b)].add(donor)
obj=collections.defaultdict(set);rows=[];prefixbad=[]
with gzip.open(P/'metadata_audit/cell_donor_ids.tsv.gz','rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  donor=x['DonorID'];prefix=donor+'_'+donor+'_'
  if not x['cell_id'].startswith(prefix):prefixbad.append(x);continue
  b=x['cell_id'][len(prefix):];obj[donor].add(b)
  if b not in raw[donor]:
   matches=sorted(index[b]);corematches=sorted(core[re.sub(r'-\d+$','',b)])
   rows.append({'DonorID':donor,'GSM':by[donor]['GSM'],'metadata_cell_id':x['cell_id'],'barcode_exact':b,'other_exact_donors':';'.join(matches),'other_exact_GSMs':';'.join(by[d]['GSM'] for d in matches),'other_exact_count':len(matches),'same_donor_after_suffix_removal':donor in corematches,'all_core_donors':';'.join(corematches),'classification':'absent_all_48_raw_sets' if not matches else ('one_other_raw_set_barcode_only' if len(matches)==1 else 'multiple_other_raw_sets_barcode_only'),'identity_reassigned':False})
def write(name,rs):
 with (O/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
write('unmatched_4407_crosssample.csv',rows);write('raw_barcode_member_versions.csv',members)
overlap=[]
for d,ss in obj.items():
 for e,tt in raw.items():overlap.append({'object_donor':d,'raw_donor':e,'raw_GSM':by[e]['GSM'],'object_cells':len(ss),'raw_cells':len(tt),'exact_shared_barcodes':len(ss&tt),'is_named_pair':d==e})
write('all_48x48_barcode_overlaps.csv',overlap)
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'object_cells':sum(map(len,obj.values())),'invalid_double_donor_prefix_count':len(prefixbad),'unmatched_count':len(rows),'class_counts':dict(collections.Counter(x['classification'] for x in rows)),'suffix_removal_rescues_same_donor':sum(x['same_donor_after_suffix_removal'] for x in rows),'raw_suffixes':sorted(set(x['suffixes'] for x in members)),'object_barcode_pattern_counts':dict(collections.Counter('16ACGT-1' if re.fullmatch('[ACGT]{16}-1',b) else 'other' for ss in obj.values() for b in ss)),'no_identity_reassignment':True,'caution':'10x nucleotide barcodes recur across libraries; any cross-sample overlap is a candidate search result, not donor identity evidence.'}
(O/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
