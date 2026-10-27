from pathlib import Path
import json,hashlib,numpy as np,pandas as pd,pyarrow.feather as feather
from pyscenic.utils import load_motif_annotations
P=Path(__file__).parent;O=P/'F5_GSE244832_regulon_v1';R=P.parent/'reference_cisTarget_hg38_mc9nr_v3';A=P.parent/'reference_cisTarget_aux_v1';genes=set(pd.read_csv(O/'filtered_genes.tsv',sep='\t').gene);rows=[];motifsets=[]

# Validate the GRN output before reference parsing and motif pruning.
adj=pd.read_csv(O/'adjacencies.tsv',sep='\t');assert set(['TF','target','importance']).issubset(adj.columns)
tfset=set((P/'F5_regulatory_resource_audit_v1/allTFs_hg38.txt').read_text().split())
assert len(adj)>0 and adj[['TF','target','importance']].notna().all().all()
assert set(adj.TF)<=tfset & genes and set(adj.target)<=genes
assert not adj[['TF','target']].duplicated().any() and (adj.TF!=adj.target).all()
assert np.isfinite(adj.importance).all() and (adj.importance>0).all()
(O/'GRN_output_validation.json').write_text(json.dumps({'status':'passed','edges':len(adj),'TFs':adj.TF.nunique(),'targets':adj.target.nunique(),'filtered_expression_genes':len(genes),'targets_without_inferred_edges':sorted(genes-set(adj.target)),'TF_and_target_ID_alignment':True,'finite_positive_importance':True,'no_duplicate_or_self_edges':True,'adjacency_sha256':hashlib.sha256((O/'adjacencies.tsv').read_bytes()).hexdigest()},indent=2))

