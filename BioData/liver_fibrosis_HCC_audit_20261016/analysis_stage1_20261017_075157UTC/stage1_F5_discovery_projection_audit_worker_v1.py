from pathlib import Path
import json,hashlib,pandas as pd,numpy as np
P=Path(__file__).parent;I=P/'F5_GSE244832_regulon_v1';O=P/'F5_discovery_projection_audit_v1';O.mkdir(exist_ok=False)
def sha(q):
 h=hashlib.sha256()
 with open(q,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
frozen=json.loads((I/'input_validation.json').read_text());hashes={}
for n in ['discovery_raw_counts.tsv','projection_raw_counts.tsv','cell_metadata.tsv','analysis_rules_frozen.json']:
 hashes[n]=sha(I/n);assert hashes[n]==frozen['files'][n]['sha256']
m=pd.read_csv(I/'cell_metadata.tsv',sep='\t',index_col=0,keep_default_na=False);x=pd.read_csv(I/'discovery_raw_counts.tsv',sep='\t',index_col=0);y=pd.read_csv(I/'projection_raw_counts.tsv',sep='\t',index_col=0);a=pd.read_csv(I/'aucell_matrix.csv',index_col=0);d=m.regulon_role.eq('discovery_MASH');assert x.index.equals(m.index[d]) and y.index.equals(m.index) and a.index.equals(y.index);assert x.columns.equals(y.columns);assert np.array_equal(x.to_numpy(),y.loc[x.index].to_numpy());assert len(set(x.index)-set(y.index))==0 and len(x)==2022 and len(y)==3733
assert m.loc[d,'condition_Excel_original'].eq('MASH').all();assert not m.donor.isin(['D9','D10']).any();assert set(m.subcluster_preliminary.astype(str))<=set(map(str,range(8)))
tables={}
for name,frame in [('discovery',m[d]),('projection_including_discovery',m),('additional_projection_only',m[~d])]:
 z=frame.groupby(['condition_Excel_original','F_stage_from_Figure4C'],dropna=False).agg(cells=('donor','size'),donors=('donor','nunique'));z.to_csv(O/(name+'_strata.tsv'),sep='\t');tables[name]={'cells':len(frame),'donors':frame.donor.nunique(),'donor_ids':sorted(frame.donor.unique()),'strata':z.reset_index().to_dict('records')}
m.groupby(['regulon_role','donor','condition_Excel_original','F_stage_from_Figure4C','subcluster_preliminary']).size().rename('cells').to_csv(O/'donor_cluster_design.tsv',sep='\t');m.assign(in_discovery=d).to_csv(O/'cell_ID_membership.tsv.gz',sep='\t')
commands={}
for stage in ['GRNBoost2','motif_pruning','AUCell']:
 c=json.loads((I/(stage+'_complete.json')).read_text())['command'];commands[stage]=c
 assert str(I/('projection_raw_counts.tsv' if stage=='AUCell' else 'discovery_raw_counts.tsv')) in c
v={'status':'passed','sets':tables,'same_8189_genes_in_both_inputs':len(x.columns)==8189,'discovery_IDs_exact_subset_of_projection':True,'raw_count_values_identical_on_discovery_subset':True,'AUCell_ID_order_matches_predeclared_projection':True,'frozen_input_hashes_unchanged':hashes,'frozen_rules_utc':json.loads((I/'analysis_rules_frozen.json').read_text())['utc'],'primary_CSI_and_RSS_cells':2022,'primary_CSI_and_RSS_donors':8,'identity':'Reviewed macrophage-like numeric subclusters0-7;not all myeloid cells;unresolved8 excluded','interpretation':'Predeclared same-cohort projection; not an expanded GRN discovery input, not independent external validation','actual_completed_commands':commands};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
