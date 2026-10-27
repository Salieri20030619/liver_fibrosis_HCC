from pathlib import Path
import urllib.request,json,hashlib,datetime,traceback,collections
import pandas as pd,numpy as np
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'F3F4_KEGG_function_v2';O.mkdir(exist_ok=False);op=urllib.request.build_opener(urllib.request.ProxyHandler({}));audit=[]
def bh(p):
 p=np.array(p);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
try:
 for n in ['genes.tsv','links.tsv','pathways.tsv']:
  b=(P/'KEGG_official_sources_v2'/n).read_bytes();(O/n).write_bytes(b)
 sym=collections.defaultdict(set)
 for line in (O/'genes.tsv').read_text().splitlines():
  a=line.split('\t');symbol=a[-1].split(';')[0].split(',')[0].strip()
  if symbol:sym[symbol].add(a[0])
 by={k:next(iter(v)) for k,v in sym.items() if len(v)==1};paths=collections.defaultdict(set)
 for line in (O/'links.tsv').read_text().splitlines():
  gene,path=line.split('\t');paths[path.replace('path:','')].add(gene)
 names=dict(line.split('\t',1) for line in (O/'pathways.tsv').read_text().splitlines());results=[]
 for f in (P/'F3F4_donor_analysis_v1').glob('*/all_gene_donor_rank_sum.tsv.gz'):
  x=pd.read_csv(f,sep='\t');x['KEGG_ID']=x.gene_symbol.map(by);x=x[x.KEGG_ID.notna()];duplicated=x.KEGG_ID.duplicated(keep=False);amb=int(duplicated.sum());x=x[~duplicated];universe=set(x.KEGG_ID);d=O/f.parent.name;d.mkdir();x[['gene_id','gene_symbol','KEGG_ID']].to_csv(d/'tested_background_mapping.tsv',sep='\t',index=False)
  for direction in ['up','down']:
   sig=x[(x.p_value<.05)&((x.mean_donor_log2CPM_difference>0) if direction=='up' else (x.mean_donor_log2CPM_difference<0))];selected=set(sig.KEGG_ID);rows=[]
   for term,genes in paths.items():
    bg=genes&universe
    if not bg:continue
    hit=bg&selected;rows.append({'pathway':term,'name':names.get(term,term),'tested_background':len(universe),'pathway_background':len(bg),'selected_genes':len(selected),'overlap':len(hit),'p_value':float(hypergeom.sf(len(hit)-1,len(universe),len(bg),len(selected))) if selected else 1.,'overlap_gene_IDs':';'.join(sorted(hit))})
   t=pd.DataFrame(rows);t['BH_q']=bh(t.p_value);t=t.sort_values('p_value');t.to_csv(d/(direction+'_all_KEGG_ORA.tsv'),sep='\t',index=False);show=t.head(15).iloc[::-1];fig,ax=plt.subplots(figsize=(10,7));ax.barh(show.name.str.replace(' - Homo sapiens (human)','',regex=False),-np.log10(np.maximum(show.p_value,1e-300)));ax.set_xlabel('-log10 pathway p (full table includes BH q)');ax.set_title(f.parent.name+' | '+direction+' | exploratory nominal-DE ORA');fig.tight_layout();fig.savefig(d/(direction+'_F3F_F4_function.png'),dpi=160);plt.close(fig);results.append({'comparison':f.parent.name,'direction':direction,'mapped_universe':len(universe),'ambiguous_duplicate_ID_rows_excluded':amb,'nominal_selected':len(selected),'tested_pathways':len(t),'BH_q_below_05':int((t.BH_q<.05).sum()),'top3':t.head(3).to_dict('records')})
 (O/'analysis_design.json').write_text(json.dumps({'panels':'F3F and functional component of F4D-G; not full GSVA/metabolic panel reproduction','method':'KEGG over-representation, hypergeometric upper tail, up/down separately, raw gene p<0.05 from baseline significance convention; BH across all tested pathways per direction. No fold-change or gene-set size cutoff. Background all tested genes mapped by exact unique first KEGG symbol. Ambiguous repeated KEGG IDs omitted and counted.','interpretation':'Exploratory only: no individual gene passes BH q<0.05 in current donor analyses. Pathway enrichment of nominal gene lists does not establish mechanism, lineage or axis. No arbitrary pathway selection based on expected axis.','source':'Official KEGG REST; current dated snapshot, not known original article version'},indent=2));(O/'results.json').write_text(json.dumps(results,indent=2));print(json.dumps({'status':'complete','contrasts_directions':len(results)}),flush=True)
except Exception:(O/'error.txt').write_text(traceback.format_exc());raise
