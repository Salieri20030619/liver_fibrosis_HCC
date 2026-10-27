import pathlib,tarfile,io,json,gzip,hashlib,datetime
import h5py,numpy as np,pandas as pd
from PIL import Image
from scipy.io import mmread
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016');p=r/'GSM7661256_HCC2R.tar.gz'
with tarfile.open(p) as t:
 def b(n):return t.extractfile('HCC2R/'+n).read()
 with h5py.File(io.BytesIO(b('filtered_feature_bc_matrix.h5'))) as h:
  m=h['matrix'];shape=m['shape'][:].tolist();bc=[x.decode() for x in m['barcodes'][:]];genes=[x.decode() for x in m['features/id'][:]];data=m['data'][:];ind=m['indices'][:];ptr=m['indptr'][:]
  genome=sorted(set(x.decode() for x in m['features/genome'][:]))
 pos=pd.read_csv(io.BytesIO(b('spatial/tissue_positions_list.csv')),header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']);inside=pos[pos.in_tissue==1]
 mtx=mmread(gzip.GzipFile(fileobj=io.BytesIO(b('filtered_feature_bc_matrix/matrix.mtx.gz')))).tocsc();bc2=gzip.decompress(b('filtered_feature_bc_matrix/barcodes.tsv.gz')).decode().splitlines();genes2=[x.split('\t')[0] for x in gzip.decompress(b('filtered_feature_bc_matrix/features.tsv.gz')).decode().splitlines()]
 images={}
 for name in ['tissue_lowres_image.png','tissue_hires_image.png']:
  raw=b('spatial/'+name);im=Image.open(io.BytesIO(raw));im.load();images[name]={'size':im.size,'sha256':hashlib.sha256(raw).hexdigest()}
  if 'lowres' in name:(r/'GSM7661256_HCC2R_lowres.png').write_bytes(raw)
 v={'resource':'GSE238264','sample':'GSM7661256/HCC2R','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'file':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'shape':shape,'nnz':len(data),'genome':genome,'duplicate_genes':len(genes)-len(set(genes)),'duplicate_barcodes':len(bc)-len(set(bc)),'negative':int((data<0).sum()),'nonfinite':int((~np.isfinite(data)).sum()),'positions':len(pos),'in_tissue_positions':len(inside),'matrix_without_coordinates':len(set(bc)-set(pos.barcode)),'in_tissue_without_matrix':len(set(inside.barcode)-set(bc)),'duplicate_position_barcodes':int(pos.barcode.duplicated().sum()),'h5_mtx_same_shape':list(mtx.shape)==shape,'h5_mtx_same_barcodes':bc==bc2,'h5_mtx_same_gene_order':genes==genes2,'h5_mtx_same_sparse_data':bool(np.array_equal(data,mtx.data) and np.array_equal(ind,mtx.indices) and np.array_equal(ptr,mtx.indptr)),'images':images,'scalefactors':json.loads(b('spatial/scalefactors_json.json')),'scope':'One treated patient only; six other GEO samples not downloaded'}
 (r/'GSE238264_HCC2R_full_validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
