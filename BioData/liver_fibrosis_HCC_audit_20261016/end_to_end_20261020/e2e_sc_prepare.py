from pathlib import Path
import tarfile,gzip,io,json,traceback,subprocess,os
import pandas as pd,numpy as np,scipy.sparse as sp
from scipy.io import mmwrite
import anndata as ad
P=Path(__file__).parent;O=P/'GSE202379_core39_QC';O.mkdir(exist_ok=False)
try:
 meta=pd.read_csv(P/'GSE202379_core39.tsv',sep='\t').set_index('sample_id');assert len(meta)==39 and meta.donor_label.is_unique
 blocks=[];obs=[];summ=[];genes=None
 with tarfile.open(P.parent/'freeze_priority_20261019/GSE202379_RAW.tar') as tar:
  for mem in tar:
   gsm=mem.name.split('_')[0]
   if gsm not in meta.index or 'raw_counts' not in mem.name:continue
   x=pd.read_csv(io.BytesIO(gzip.decompress(tar.extractfile(mem).read())),index_col=0)
   if genes is None:genes=x.index.astype(str)
   assert genes.equals(x.index.astype(str)) and x.columns.is_unique
   a=sp.csr_matrix(x.to_numpy(dtype=np.float32).T);total=np.asarray(a.sum(axis=1)).ravel();ng=np.asarray((a>0).sum(axis=1)).ravel();mt=genes.str.upper().str.startswith('MT-');pct=np.asarray(a[:,mt].sum(axis=1)).ravel()/np.maximum(total,1)*100
   metrics={'UMI':total,'genes':ng,'percent_mt':pct};keep=total>0;limits={}
   for name,v in metrics.items():
    med=float(np.median(v));mad=float(np.median(np.abs(v-med)));limits[name]=[med-3*mad,med+3*mad]
    if mad>0:keep&=(v>=med-3*mad)&(v<=med+3*mad)
   cellids=[gsm+'|'+str(c) for c in x.columns[keep]];m=pd.DataFrame({'donor':meta.loc[gsm,'donor_label'],'GSM':gsm,'F':int(meta.loc[gsm,'F_direct']),'study':'GSE202379','UMI':total[keep],'genes':ng[keep],'percent_mt':pct[keep]},index=cellids)
   blocks.append(a[keep]);obs.append(m);summ.append({'GSM':gsm,'donor':meta.loc[gsm,'donor_label'],'F':int(meta.loc[gsm,'F_direct']),'raw_cells':len(x.columns),'QC_pass':int(keep.sum()),'mt_features':int(mt.sum()),'thresholds':json.dumps(limits)})
   print('QC',gsm,a.shape,int(keep.sum()),flush=True)
 assert len(blocks)==39
 a=ad.AnnData(sp.vstack(blocks).tocsr(),obs=pd.concat(obs),var=pd.DataFrame(index=genes));assert a.obs_names.is_unique
 a.write_h5ad(O/'core39_counts.h5ad',compression='gzip');pd.DataFrame(summ).to_csv(O/'per_donor_QC.tsv',sep='\t',index=False)
 mmwrite(O/'counts.mtx',a.X.T);a.obs.to_csv(O/'cell_metadata.tsv',sep='\t');pd.Series(genes).to_csv(O/'genes.tsv',sep='\t',index=False,header=False)
 v={'status':'completed','shape':list(a.shape),'donors':39,'F_counts':meta.F_direct.value_counts().to_dict(),'raw_cells':sum(v['raw_cells'] for v in summ),'method':'Original 3MAD on UMI/genes/MT per source sample, symmetric bounds. Zero MAD does not remove cells. Threshold details reconstructed; doublet removal not specified in paper.','mt_features':int(mt.sum())};(O/'validation.json').write_text(json.dumps(v,indent=2));print('COMPLETE',json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
