from pathlib import Path
import pandas as pd,numpy as np,json,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;src=P/'F2_Roe_F3_F4_subclusters_v1';out=P/'F3_F4_review_plots_v1';out.mkdir(exist_ok=False);results={}
for d in src.iterdir():
 if not d.is_dir() or not (d/'result.json').exists():continue
 c=pd.read_csv(d/'subcluster_coordinates_metadata.tsv.gz',sep='\t',index_col=0);e=pd.read_csv(d/'subcluster_marker_evidence.tsv',sep='\t');fig,ax=plt.subplots(figsize=(8,6))
 for cl,q in c.groupby('subcluster_preliminary'):ax.scatter(q.UMAP1,q.UMAP2,s=2,label=str(cl),rasterized=True)
 ax.legend(title='Preliminary cluster',bbox_to_anchor=(1,1),markerscale=4);ax.set_title(d.name+' | identities pending review');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');fig.tight_layout();fig.savefig(out/(d.name+'_UMAP.png'),dpi=180);plt.close(fig)
 genes=list(dict.fromkeys(e.marker));clusters=sorted(e.subcluster.unique());fig,ax=plt.subplots(figsize=(14,6));xx=[];yy=[];ss=[];cc=[]
 for _,r in e.iterrows():xx.append(genes.index(r.marker));yy.append(clusters.index(r.subcluster));ss.append(100*r.fraction_detected);cc.append(r.mean_log1p_cp10k)
 im=ax.scatter(xx,yy,s=ss,c=cc,cmap='viridis');ax.set_xticks(range(len(genes)),genes,rotation=90);ax.set_yticks(range(len(clusters)),clusters);ax.set_title(d.name+' | area=detection fraction; color=mean log1p CP10k');fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(out/(d.name+'_marker_dotplot.png'),dpi=180);plt.close(fig);results[d.name]={'cells':len(c),'clusters':len(clusters),'fine_identities_assigned':False}
(out/'result.json').write_text(json.dumps(results,indent=2));print(json.dumps(results))
