from pathlib import Path
import csv,json,gzip,tarfile,collections,time
import numpy as np
from scipy.io import mmread
from threadpoolctl import threadpool_limits
P=Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/zenodo_seurat_member_20261020');O=P/'barcode_followup'
req=list(csv.DictReader((O/'count_comparison_requests.csv').open()));cells=(O/'object_count_cells.txt').read_text().splitlines();genes=(O/'object_count_features.txt').read_text().splitlines();ci={x:i for i,x in enumerate(cells)}
with threadpool_limits(limits=1):obj=mmread(O/'object_selected_counts.mtx').tocsc()
rows=[];mappingrows=[]
with tarfile.open('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016/freeze_priority_20261019/GSE289173_RAW.tar','r:') as t:
 for donor in sorted(set(x['raw_donor'] for x in req)):
  rr=[x for x in req if x['raw_donor']==donor];gsm=rr[0]['raw_GSM'];prefix=gsm+'_'+donor
  f=[x.split('\t') for x in gzip.decompress(t.extractfile(prefix+'_features.tsv.gz').read()).decode().splitlines()];barcodes=gzip.decompress(t.extractfile(prefix+'_barcodes.tsv.gz').read()).decode().splitlines();bi={x:i for i,x in enumerate(barcodes)}
  ids=[x[0] for x in f];symbols=[x[1] for x in f];idcounts=collections.Counter(ids);symcounts=collections.Counter(symbols)
  idmap={x:i for i,x in enumerate(ids) if idcounts[x]==1};symmap={x:i for i,x in enumerate(symbols) if symcounts[x]==1}
  idhits=sum(x in idmap for x in genes);symhits=sum(x in symmap for x in genes);mode='exact_unique_gene_id' if idhits>=symhits else 'exact_unique_gene_symbol';gm=idmap if idhits>=symhits else symmap
  oi=[i for i,x in enumerate(genes) if x in gm];ri=[gm[genes[i]] for i in oi]
  mappingrows.append({'raw_donor':donor,'raw_GSM':gsm,'mode':mode,'object_features':len(genes),'raw_features':len(f),'aligned_unambiguous_features':len(oi),'object_features_unaligned':len(genes)-len(oi),'raw_duplicate_symbols_excluded':sum(c for c in symcounts.values() if c>1)})
  with gzip.GzipFile(fileobj=t.extractfile(prefix+'_matrix.mtx.gz')) as stream:
   with threadpool_limits(limits=1):raw=mmread(stream).tocsc()
  for r in rr:
   a=obj[oi,ci[r['object_cell_id']]].tocsc();b=raw[ri,bi[r['barcode']]].tocsc();delta=a-b;delta.eliminate_zeros()
   rows.append({**r,'gene_alignment_mode':mode,'aligned_genes':len(oi),'object_nonzero_aligned':a.nnz,'raw_nonzero_aligned':b.nnz,'object_sum_aligned':float(a.sum()),'raw_sum_aligned':float(b.sum()),'differing_gene_counts':delta.nnz,'absolute_count_difference':float(abs(delta).sum()),'exact_equal_on_aligned_genes':len(oi)>1000 and a.nnz>0 and delta.nnz==0,'identity_reassigned':False})
  print(donor,len(rr),'checked',flush=True)
for name,rs in [('count_vector_comparison.csv',rows),('gene_alignment.csv',mappingrows)]:
 with (O/name).open('w',newline='') as out:
  w=csv.DictWriter(out,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
summary={'pairs_compared':len(rows),'control_pairs':sum(x['type']=='named_pair_control' for x in rows),'control_exact_equal':sum(x['type']=='named_pair_control' and x['exact_equal_on_aligned_genes'] for x in rows),'cross_pairs':sum(x['type']=='cross_sample_candidate' for x in rows),'cross_exact_equal':sum(x['type']=='cross_sample_candidate' and x['exact_equal_on_aligned_genes'] for x in rows),'min_aligned_genes':min(x['aligned_genes'] for x in rows),'max_aligned_genes':max(x['aligned_genes'] for x in rows),'no_identity_reassignment':True}
(O/'count_comparison_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
