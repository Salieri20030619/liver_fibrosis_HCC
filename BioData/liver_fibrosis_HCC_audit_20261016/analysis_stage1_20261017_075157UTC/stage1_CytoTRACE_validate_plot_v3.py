from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;O=P/'F3_CytoTRACE_validated_v3';O.mkdir(exist_ok=False);results=[]
for g in ['GSE136103','GSE244832']:
 src=P/'F3_CytoTRACE_official_v2'/g;d=P/'F3_DDRTree_unrooted_v1'/g;out=O/g;out.mkdir();scores=pd.read_csv(src/'cell_scores.tsv',sep='\t',index_col=0);meta=pd.read_csv(d/'cell_metadata.tsv',sep='\t',index_col=0);xy=pd.read_csv(d/'DDRTree_coordinates.tsv',sep='\t',index_col=0);assert scores.index.equals(meta.index) and xy.index.equals(meta.index);ok=np.isfinite(scores.CytoTRACE);assert scores.CytoTRACE[ok].between(0,1).all();full=meta.join(scores);full['group']=np.where(full.F_stage_from_Figure4C.astype(str).isin(['3','4','3.0','4.0']),'F3_F4','F0') if g=='GSE244832' else full['disease status'];eligible=full.disease_comparison_eligible.astype(str).str.lower().eq('true');don=full[eligible & ok].groupby(['donor','sorting','group'],observed=True).CytoTRACE.agg(['mean','count']).reset_index();don.to_csv(out/'donor_scores.tsv',sep='\t',index=False)
 fig,axs=plt.subplots(1,2,figsize=(13,5));im=axs[0].scatter(xy.iloc[:,0][ok],xy.iloc[:,1][ok],c=scores.CytoTRACE[ok],s=5,cmap='viridis',vmin=0,vmax=1);fig.colorbar(im,ax=axs[0],label='Official CytoTRACE');axs[0].set(xlabel='DDRTree 1',ylabel='DDRTree 2',title='Score overlay; no root or direction assigned');keys=[];vals=[]
 for (sort,group),z in don.groupby(['sorting','group']):keys.append(str(sort)+'\n'+group);vals.append(z['mean'].to_numpy())
 axs[1].boxplot(vals,tick_labels=keys)
 for i,v in enumerate(vals):axs[1].scatter(np.full(len(v),i+1),v,s=20)
 axs[1].set(ylabel='Mean CytoTRACE per donor / sorting',title='Descriptive within-cohort scores');fig.suptitle(g+' | Broad mesenchymal compartment, not confirmed HSC lineage');fig.tight_layout();fig.savefig(out/'F3G_CytoTRACE_with_unrooted_topology.png',dpi=160);plt.close(fig)
 r={'cohort':g,'input_cells':len(meta),'scored_finite_cells':int(ok.sum()),'unscored_cells':int((~ok).sum()),'ids_exact':True,'score_range':[float(scores.CytoTRACE.min()),float(scores.CytoTRACE.max())],'score_vs_detected_genes_spearman':float(scores.CytoTRACE.corr(scores.detected_genes,method='spearman')),'donor_sorting_rows':len(don),'independent_donors':don.donor.nunique(),'root_assigned':False,'interpretation':'Within-cohort differentiation-potential estimate; no HCC mechanism or lineage causality. Gene-detection and sorting effects require review.'};results.append(r);(out/'validation.json').write_text(json.dumps(r,indent=2))
(O/'validation.json').write_text(json.dumps(results,indent=2));print(json.dumps(results))
