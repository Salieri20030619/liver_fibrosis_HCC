from pathlib import Path
import json,pandas as pd,pyarrow as pa
P=Path(__file__).parent;I=P/'F5_GSE244832_regulon_v1';R=P.parent/'reference_cisTarget_hg38_mc9nr_v3';g=pd.read_csv(I/'filtered_genes.tsv',sep='\t');out=g.copy();counts={}
for n,q in enumerate(sorted(R.glob('*.feather'))):
 with pa.memory_map(str(q),'r') as src:genes=set(pa.ipc.open_file(src).schema.names)-{'motifs'}
 col='present_in_'+('10kb_database' if '10kb' in q.name else 'proximal_database');out[col]=out.gene.isin(genes);counts[q.name]={'overlap':int(out[col].sum()),'missing':int((~out[col]).sum())}
tf=set((P/'F5_regulatory_resource_audit_v1/allTFs_hg38.txt').read_text().split());out['listed_TF']=out.gene.isin(tf);out.to_csv(I/'reference_gene_coverage.tsv',sep='\t',index=False);v={'status':'verified','discovery_genes':len(g),'database_coverage':counts,'gene_alias_remapping':False,'interpretation':'Missing database genes are explicitly listed; motif support unavailable for those genes. This is an older official hg38/refseq-r80 mc9nr reference, not a newly re-annotated genome.'};(I/'reference_gene_coverage.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
