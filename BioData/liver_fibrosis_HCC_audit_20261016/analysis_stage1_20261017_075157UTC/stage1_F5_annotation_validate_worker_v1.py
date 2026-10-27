from pathlib import Path
import json,hashlib
from pyscenic.utils import load_motif_annotations
P=Path(__file__).parent;A=P.parent/'reference_cisTarget_aux_v1';f=A/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl';s=json.loads((A/'complete.json').read_text());assert f.stat().st_size==s['annotation_bytes'] and hashlib.sha256(f.read_bytes()).hexdigest()==s['annotation_sha256'];a=load_motif_annotations(str(f));assert len(a)>1000
v={'status':'parsed','bytes':f.stat().st_size,'sha256':s['annotation_sha256'],'filtered_annotation_rows':len(a),'index_names':list(a.index.names),'unique_each_index_level':{str(n):a.index.get_level_values(i).nunique() for i,n in enumerate(a.index.names)},'columns':list(a.columns),'official_parser':'pyscenic.utils.load_motif_annotations0.12.1; default orthologyFDRfilters'};(A/'parser_validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
