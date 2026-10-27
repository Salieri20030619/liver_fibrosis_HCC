import pathlib,json,datetime,hashlib,csv,collections,re
import h5py,numpy as np
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');p=r/'BIT_range_resume/Liver.h5ad.complete';o=r/'BIT_final_mapping_audit';o.mkdir(exist_ok=False)
v={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'url':'https://zenodo.org/api/records/14231098/files/Input_Data.tar.gz/content','file':str(p),'bytes':p.stat().st_size,'expected_member_bytes':4829857046,'archive_fully_downloaded_or_validated':False}
h=hashlib.sha256()
with p.open('rb') as f:
 for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
v['sha256']=h.hexdigest();assert v['sha256']=='a77b87999ae2ba3e08d4c04b04f30f214c0a5a404467e036a5c536790ab008e3'
dec=lambda x:x.decode() if isinstance(x,bytes) else str(x)
with h5py.File(p,'r') as f:
 obs=f['obs'];ids=[dec(x) for x in obs['index'][:]];cats=[dec(x) for x in obs['cell_type/categories'][:]];labels=[cats[i] for i in obs['cell_type/codes'][:]]
 v['obs_fields']=list(obs);v['celltype_counts']=dict(collections.Counter(labels));v['obs_duplicate_IDs']=len(ids)-len(set(ids))
 tf=list(csv.DictReader((r/'GSE227265_TF_cell_ID_order.tsv').open(),delimiter='\t'));tfids=[x['cell_id'] for x in tf];lookup={x:i for i,x in enumerate(ids)}
 v['TF_ID_alignment']={'TF_cells':len(tfids),'BIT_cells':len(ids),'intersection':len(set(ids)&set(tfids)),'same_set':set(ids)==set(tfids),'same_order':ids==tfids}
 with (o/'TF_to_BIT_exact_ID_join.tsv').open('x',newline='') as out:
  w=csv.writer(out,delimiter='\t');w.writerow(['TF_matrix_column_1based','TF_physical_CSV_column_1based','cell_id','BIT_obs_row_1based','provided_cell_type'])
  for j,cell in enumerate(tfids):i=lookup[cell];w.writerow([j+1,j+2,cell,i+1,labels[i]])
 v['matrices']={}
 for name in ['X','obsm/fragment_paired']:
  g=f[name];shape=list(map(int,g.attrs['shape']));ip=g['indptr'][:];n=len(g['data']);badidx=nonfinite=negative=zeros=0
  for st in range(0,n,4000000):
   ix=g['indices'][st:st+4000000];dat=g['data'][st:st+4000000];badidx+=int(((ix<0)|(ix>=shape[1])).sum());nonfinite+=int((~np.isfinite(dat)).sum());negative+=int((dat<0).sum());zeros+=int((dat==0).sum())
  v['matrices'][name]={'shape':shape,'stored_entries':n,'indices_length':len(g['indices']),'indptr_length':len(ip),'indptr_start':int(ip[0]),'indptr_end':int(ip[-1]),'indptr_nondecreasing':bool(np.all(np.diff(ip)>=0)),'index_out_of_range':badidx,'nonfinite':nonfinite,'negative':negative,'stored_zero_entries':zeros,'empty_rows':int((np.diff(ip)==0).sum()),'full_values_and_indices_scanned':True}
  v['matrices'][name]['structural_pass']=len(ip)==shape[0]+1 and ip[0]==0 and ip[-1]==n==len(g['indices']) and bool(np.all(np.diff(ip)>=0)) and badidx==nonfinite==negative==0
 ref=f['uns/reference_sequences'];refs={dec(k):int(n) for k,n in zip(ref['reference_seq_name'][:],ref['reference_seq_length'][:])};var=f['var'];seen=collections.Counter();widths=collections.Counter();bad=0;selected=int(var['selected'][:].sum());prev={};gaps=0
 for st in range(0,len(var['_index']),100000):
  for raw in var['_index'][st:st+100000]:
   s=dec(raw);m=re.fullmatch(r'([^:]+):(\d+)-(\d+)',s)
   if not m:bad+=1;continue
   chrom,a,b=m.group(1),int(m.group(2)),int(m.group(3));seen[chrom]+=1;widths[b-a]+=1
   if chrom not in refs or a<0 or b<=a or b>refs[chrom]:bad+=1
   if a!=prev.get(chrom,0):gaps+=1
   prev[chrom]=b
 v['var']={'rows':len(var['_index']),'fields':list(var),'selected_true':selected,'invalid_coordinates_vs_embedded_reference':bad,'width_counts':dict(widths),'chromosome_counts':dict(seen),'discontinuities':gaps,'reference_lengths':refs,'last_coordinate_end_per_chromosome':prev,'coordinate_reference':'Embedded reference names and lengths; no explicit assembly label asserted'}
 v['other_annotation_keys']={name:list(f[name]) for name in ['uns','obsm','obsp','varm','varp']}
 v['mapping_conclusion']={'celltype_to_TF':'Exact ID join recovered for all 18631 TF columns; reorder explicitly. Supplied labels not independently biologically validated.','patient_donor_HCC_CCA':'Not recovered; no explicit donor/patient/GSM or HCC/CCA assignment in obs or other annotation containers. Barcode suffixes not used as donor evidence.','original_GEO_peak_rows':'Not recovered: original GEO MTX has 1089114 peaks, this object has 6062095 genomic bins. Neither row ordering nor equivalence may be inferred from dimensions or selected mask.','decision':'Conditional inclusion for celltype-labelled TF exploratory summaries; not eligible for donor-aware HCC-only inference or original peak-coordinate analyses until missing mappings are resolved.'}
v['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
with (o/'BIT_integrity_mapping_final.json').open('x') as out:json.dump(v,out,indent=2)
print(json.dumps(v),flush=True)
