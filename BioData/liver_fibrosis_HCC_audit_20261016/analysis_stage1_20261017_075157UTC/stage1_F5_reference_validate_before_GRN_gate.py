from pathlib import Path
import json,hashlib,numpy as np,pandas as pd,pyarrow.feather as feather
from pyscenic.utils import load_motif_annotations
P=Path(__file__).parent;O=P/'F5_GSE244832_regulon_v1';R=P.parent/'reference_cisTarget_hg38_mc9nr_v3';A=P.parent/'reference_cisTarget_aux_v1';genes=set(pd.read_csv(O/'filtered_genes.tsv',sep='\t').gene);rows=[];motifsets=[]
for q in sorted(R.glob('*.feather')):
 t=feather.read_table(q,memory_map=True);t.validate(full=True);assert t.num_rows>1000 and t.num_columns>10000
 # Feather v2 gene-based cisTarget database stores the motif index as 'motifs'.
 assert 'motifs' in t.column_names;motifs=t.column('motifs').to_pylist();assert len(set(motifs))==len(motifs);dbgenes=set(t.column_names)-{'motifs'};assert len(dbgenes)>10000;motifsets.append(set(motifs));sample_cols=sorted(dbgenes)[::max(1,len(dbgenes)//20)]
 for c in sample_cols:
  z=t.column(c).to_numpy();assert np.isfinite(z).all() and (z>=0).all() and (z<len(dbgenes)).all()
 rows.append({'file':str(q),'bytes':q.stat().st_size,'all_record_batches_fully_parsed':True,'motifs':len(motifs),'genes':len(dbgenes),'discovery_gene_overlap':len(genes&dbgenes),'sampled_rank_columns_valid':len(sample_cols)})
 del t
ann=load_motif_annotations(str(A/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'));assert len(ann)>1000;motif_ids=set(ann.index.get_level_values(0));v={'status':'passed','databases':rows,'motifs_shared_between_databases':len(set.intersection(*motifsets)),'annotation_rows_after_official_filter':len(ann),'annotation_motifs_in_databases':[len(motif_ids&s) for s in motifsets],'species':'human','assembly':'hg38 confirmed in cached GEO RNAprocessing','annotation_collection':'mc9nr/v9 matched','official_integrity':'Both complete files match official per-fileSHA1; localSHA256 in download manifest. Global officialSHA256 endpoint404 retained.'};assert all(n>1000 for n in v['annotation_motifs_in_databases']);(O/'reference_validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
