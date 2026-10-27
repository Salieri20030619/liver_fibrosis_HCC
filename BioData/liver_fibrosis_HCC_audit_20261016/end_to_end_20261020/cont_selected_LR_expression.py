from pathlib import Path
import time,json,subprocess,numpy as np,pandas as pd,anndata as ad
P=Path(__file__).parent;O=P/'continuation_selected_LR_expression_v1';O.mkdir(exist_ok=False);L=P/'continuation_ligand_to_GSEA_v1';I=P/'continuation_network_communication_inputs_v1/CellChat'
for i in range(1440):
 if (L/'validation.json').exists():break
 time.sleep(10)
else:raise RuntimeError('Ligand selection pending4hours')
lv=json.loads((L/'validation.json').read_text());ligand=lv.get('ligand')
if not ligand:
 (O/'validation.json').write_text(json.dumps({'status':'not_applicable_negative_no_ligand'}));raise SystemExit(0)
cc=pd.read_csv(P/'continuation_CellChat_expanded_v1/all_record_state_specific_communications.tsv.gz',sep='\t');cc=cc[cc.ligand.eq(ligand)&cc.source.str.startswith('Myeloid_')&cc.target.eq('Stromal')&cc.F.ge(3)&cc.pval.le(.05)];cc.to_csv(O/'selected_ligand_edges.tsv',sep='\t',index=False)
r="db<-readRDS('"+str(I/'CellChatDB.human.rds')+"');write.table(db$complex,'"+str(O/'official_complex_subunits.tsv')+"',sep='\\t',col.names=NA,quote=FALSE)";z=subprocess.run(['/usr/bin/Rscript','--vanilla','-e',r],capture_output=True,text=True);(O/'complex_export.log').write_text(z.stdout+z.stderr);assert z.returncode==0
complexes=pd.read_csv(O/'official_complex_subunits.tsv',sep='\t',index_col=0).fillna('');features=sorted(set(cc.receptor)|{ligand});mapping=[];genes=set()
for feature in features:
 subunits=[x for x in complexes.loc[feature].tolist() if x] if feature in complexes.index else [feature]
 for gene in subunits:genes.add(gene);mapping.append({'feature':feature,'gene':gene,'mapping':'official_CellChat_complex' if feature in complexes.index else 'single_gene'})
pd.DataFrame(mapping).to_csv(O/'selected_LR_gene_subunits.tsv',sep='\t',index=False)
m=pd.read_csv(P/'continuation_annotation_refined_v1/cell_metadata_reviewed.tsv.gz',sep='\t',index_col=0);a=ad.read_h5ad(P/'continuation_annotation_evidence_v1/core55_record_counts.h5ad');assert m.index.equals(a.obs_names);keep=~m.reviewed_broad_lineage.str.contains('mixed');m=m[keep].copy();a=a[keep].copy();m['group']=m.reviewed_broad_lineage;m.loc[m.group.eq('Myeloid'),'group']=m.loc[m.group.eq('Myeloid'),'refined_cluster'];matched=sorted(genes&set(a.var_names));total=np.asarray(a.X.sum(1)).ravel();x=a[:,matched].X.toarray();normalized=np.log1p(x/np.maximum(total[:,None],1)*1e4);pos=pd.Series(np.arange(len(m)),index=m.index);rows=[]
for (study,record,F,group),g in m.groupby(['study','donor','F','group']):
 ix=pos.loc[g.index].to_numpy()
 for j,gene in enumerate(matched):rows.append({'study':study,'source_record':record,'F':F,'group':group,'gene':gene,'cells':len(ix),'eligible10cells':len(ix)>=10,'mean_log1pCP10K':float(normalized[ix,j].mean()),'fraction_detected':float((x[ix,j]>0).mean())})
pd.DataFrame(rows).to_csv(O/'record_group_LR_expression.tsv',sep='\t',index=False);v={'status':'completed','ligand':ligand,'receptors':sorted(set(cc.receptor)),'matched_subunit_genes':matched,'missing_genes':sorted(genes-set(matched)),'rows':len(rows),'scope':'Observed raw counts normalized within cell, then source-record/group means and detection fractions. Groups under10 cells retained with eligibility flag. Complexes decomposed using the exact installed CellChatDB, not guessed by punctuation.','limits':'Expression association only; cannot resolve ambientRNA, uptake or secretion; broadStromal not pureaHSC. No cell-independent pvalues or causal claim.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
