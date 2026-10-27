from pathlib import Path
import requests,json,hashlib,traceback
import pandas as pd,numpy as np,anndata as ad
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';O=P/'S7_current_candidates';O.mkdir(exist_ok=False)
try:
 genes=pd.read_csv(P/'F8_F9_real_full_cohort_v1/SHAP_importance.tsv',sep='\t').head(3).gene.tolist();s=requests.Session();s.trust_env=False;mapping=[]
 for g in genes:
  url='https://rest.ensembl.org/homology/symbol/homo_sapiens/'+g
  r=s.get(url,params={'target_species':'mus_musculus','type':'orthologues','sequence':'none','content-type':'application/json'},timeout=40);(O/(g+'_homology.json')).write_text(r.text);r.raise_for_status();data=r.json()['data'][0]['homologies'];one=[z for z in data if z['type']=='ortholog_one2one' and z['target']['species']=='mus_musculus'];assert len(one)==1,(g,len(one));eid=one[0]['target']['id']
  z=s.get('https://rest.ensembl.org/lookup/id/'+eid,params={'content-type':'application/json'},timeout=40);(O/(g+'_mouse_lookup.json')).write_text(z.text);z.raise_for_status();symbol=z.json()['display_name'];mapping.append({'human':g,'mouse':symbol,'mouse_ensembl':eid,'orthology_type':'ortholog_one2one','homology_url':r.url,'lookup_url':z.url})
 pd.DataFrame(mapping).to_csv(O/'verified_orthologs.tsv',sep='\t',index=False)
 a=ad.read_h5ad(D/'GSE166504_S1_S7_preparation_v1/myeloid_15_30w_source_counts.h5ad');total=np.asarray(a.X.sum(axis=1)).ravel();result=a.obs[['source_unit','timepoint','CellType','FileName']].copy()
 for z in mapping:
  assert z['mouse'] in a.var_names,z;count=np.asarray(a[:,z['mouse']].X.toarray()).ravel();result[z['human']+'_mouse_log1p_CP10k']=np.log1p(count/np.maximum(total,1)*1e4)
 result.to_csv(O/'all_myeloid_candidate_expression.tsv.gz',sep='\t');cols=[c for c in result if c.endswith('CP10k')];means=result.groupby(['source_unit','timepoint','CellType'])[cols].mean();means.to_csv(O/'source_unit_candidate_means.tsv',sep='\t')
 v={'status':'completed','cells':a.n_obs,'source_units':a.obs.source_unit.nunique(),'timepoint_cells':a.obs.timepoint.value_counts().to_dict(),'orthologs':mapping,'limit':'Optional mouse validation only; processed author myeloid cells.12 source units not asserted12 independent animals; no full38-library atlas or longitudinal lineage claim; no new formal figure.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
except Exception:
 (O/'error.txt').write_text(traceback.format_exc());print(traceback.format_exc(),flush=True);(O/'validation.json').write_text(json.dumps({'status':'blocked','error_file':'error.txt','note':'No symbol-case guessing or replacement of missing orthology'}));raise
