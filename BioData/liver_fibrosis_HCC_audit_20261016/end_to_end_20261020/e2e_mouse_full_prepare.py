from pathlib import Path
import gzip,csv,json,traceback
import pandas as pd,numpy as np,scipy.sparse as sp
from scipy.io import mmwrite
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'mouse38_full_QC';O.mkdir(exist_ok=False)
try:
 I=D/'gap_download';m=pd.read_csv(I/'GSE166504_cell_metadata.20220204.tsv.gz',sep='\t');m.index=m.FileName+'_'+m.CellID;parts=m.FileName.str.extract(r'(?P<fraction>NPC|Hepatocyte)_(?P<timepoint>Chow|[0-9]+weeks)_Animal(?P<animal_number>[0-9]+)_Capture(?P<capture>[0-9]+)');parts.index=m.index;m=pd.concat([m,parts],axis=1);assert m.index.is_unique and not m.fraction.isna().any();m['donor']=m.fraction+'|'+m.timepoint+'|'+m.animal_number;m['study']='GSE166504';m['GSM']=m.FileName;m['F']=np.nan
 genes=[];rows=[];total=np.zeros(len(m));ng=np.zeros(len(m));mt=np.zeros(len(m));mt_n=0
 with gzip.open(I/'GSE166504_cell_raw_counts.20220204.txt.gz','rt') as f:
  head=next(csv.reader(f,delimiter='\t'));cols=head[1:] if len(head)==len(m)+1 else head;assert cols==m.index.tolist()
  for i,line in enumerate(f):
   gene,rest=line.rstrip('\r\n').split('\t',1);gene=gene.strip('"');v=np.fromstring(rest,sep='\t');assert len(v)==len(m) and np.isfinite(v).all() and (v>=0).all() and (v%1==0).all();total+=v;ng+=(v>0)
   if gene.upper().startswith('MT-'):mt+=v;mt_n+=1
   rows.append(sp.csr_matrix(v.astype(np.float32).reshape(1,-1)));genes.append(gene)
   if i%5000==0:print('read genes',i,flush=True)
 assert len(set(genes))==len(genes);x=sp.vstack(rows).T.tocsr();del rows;m['UMI']=total;m['genes']=ng;m['percent_mt']=mt/np.maximum(total,1)*100 if mt_n else np.nan;m['mt_QC_applied']=bool(mt_n);keep=total>0;qc=[]
 for lib,ids in m.groupby('FileName').groups.items():
  ix=m.index.get_indexer(ids);ok=np.ones(len(ix),bool)
  for field in ['UMI','genes']+(['percent_mt'] if mt_n else []):
   v=m.iloc[ix][field].to_numpy();med=np.median(v);mad=np.median(abs(v-med))
   if mad>0:ok&=(v>=med-3*mad)&(v<=med+3*mad)
  keep[ix]&=ok;qc.append({'library':lib,'input_cells':len(ix),'QC_pass':int(ok.sum())})
 m=m.loc[keep];x=x[keep];mmwrite(O/'counts.mtx',x.T);pd.Series(genes).to_csv(O/'genes.tsv',sep='\t',index=False,header=False);m.to_csv(O/'cell_metadata.tsv',sep='\t');pd.DataFrame(qc).to_csv(O/'library_QC.tsv',sep='\t',index=False)
 counts=pd.crosstab([m.donor,m.timepoint],m.CellType);counts.to_csv(O/'source_unit_author_celltype_counts.tsv',sep='\t');counts.div(counts.sum(axis=1),axis=0).to_csv(O/'source_unit_author_celltype_fractions.tsv',sep='\t')
 v={'status':'completed','input_cells':82168,'QC_cells':len(m),'genes':len(genes),'libraries':m.FileName.nunique(),'source_units':m.donor.nunique(),'mt_features':mt_n,'limit':'Source units combine capture replicates; crossfraction animal pairing not assumed. Optionalmouse, not humanF/HCC cohort.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());raise
