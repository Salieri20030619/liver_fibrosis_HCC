from pathlib import Path
import sys,json,time
P=Path(__file__).parent;sys.path.insert(0,str(P/'python_deps'))
import pandas as pd,numpy as np,anndata as ad
O=P/'marker_expression_check_20261022';O.mkdir(exist_ok=False);start=time.time()
mk=pd.read_csv(P/'original_markers_20261022.tsv',sep='\t');m=pd.read_csv(P/'continuation_annotation_refined_v1/cell_metadata_reviewed.tsv.gz',sep='\t',index_col=0)
a=ad.read_h5ad(P/'continuation_annotation_evidence_v1/core55_record_counts.h5ad',backed='r');assert a.obs_names.tolist()==m.index.tolist()
vc=pd.Series(a.var_names).value_counts();mk['mapping_status']=mk.gene.map(lambda g:'exact_unique' if vc.get(g,0)==1 else 'absent' if vc.get(g,0)==0 else 'ambiguous_duplicate')
mk.to_csv(O/'original_marker_mapping.tsv',sep='\t',index=False);genes=sorted(mk.loc[mk.mapping_status.eq('exact_unique'),'gene'].unique());idx=a.var_names.get_indexer(genes);rows=[]
for cl,q in m.groupby('refined_cluster'):
 pos=m.index.get_indexer(q.index);z=a[pos,:].X.tocsr();total=np.asarray(z.sum(1)).ravel();v=z[:,idx].copy().astype(float);det=np.asarray((v>0).mean(0)).ravel();v=v.multiply((1e4/np.maximum(total,1))[:,None]).tocsr();v.data=np.log1p(v.data);avg=np.asarray(v.mean(0)).ravel()
 assert q.reviewed_broad_lineage.nunique()==1
 rows.append(pd.DataFrame({'current_cluster':cl,'current_final_annotation':q.reviewed_broad_lineage.iloc[0],'cells':len(q),'source_records':q.donor.nunique(),'gene':genes,'mean_log1pCP10K':avg,'fraction_detected':det}));print(cl,len(q),flush=True)
e=pd.concat(rows,ignore_index=True);e.to_csv(O/'current_cluster_marker_expression.tsv.gz',sep='\t',index=False)
joined=mk.merge(e,on='gene',how='left',validate='many_to_many');joined.to_csv(O/'original_to_current_marker_evidence.tsv.gz',sep='\t',index=False)
coverage=mk.groupby(['table','cluster','mapping_status']).size().rename('marker_rows').reset_index();coverage.to_csv(O/'marker_coverage.tsv',sep='\t',index=False)
a.file.close();v={'public_marker_rows':len(mk),'unique_public_genes':mk.gene.nunique(),'mapping':mk.mapping_status.value_counts().to_dict(),'unique_genes_evaluated':len(genes),'current_subclusters':e.current_cluster.nunique(),'expression_rows':len(e),'seconds':time.time()-start,'method':'Exact symbol only, duplicates/absent retained; mean log1p(CP10K) and detection on full raw gene denominator per cell. Gene-level annotation evidence, no aggregate signature score, no forced DTNA/aHSC relabel, no gene significance tests or ML.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
