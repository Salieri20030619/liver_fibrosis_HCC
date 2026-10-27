from pathlib import Path
import time,os,json,subprocess,pandas as pd
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';I=P/'continuation_network_communication_inputs_v1';O=P/'continuation_CellChat_expanded_v1';O.mkdir(exist_ok=False)
for i in range(240):
 if (I/'validation.json').exists() and (P/'continuation_annotation_review_acceptance.json').exists():break
 time.sleep(10)
else:raise RuntimeError('New annotation not reviewed or inputs incomplete')
r=(D/'F6_CellChat_GSE244832_v1/run_CellChat.R').read_text();r=r.replace("o<-'"+str(D/'F6_CellChat_GSE244832_v1')+"'","o<-'"+str(I/'CellChat')+"'")
# Inputs contain state-specific reviewed labels. Reuse installed original CellChat method unchanged.
r+='\nwriteLines("completed","'+str(O/'R_COMPLETE.txt')+'")\n';(O/'run.R').write_text(r);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
with open(O/'R.log','w') as f:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(O/'run.R')],env=env,stdout=f,stderr=subprocess.STDOUT)
audit=pd.read_csv(I/'CellChat/record_eligibility.tsv',sep='\t');rows=[];tabs=[]
for q in sorted((I/'CellChat').iterdir()):
 if not q.is_dir():continue
 row={'source_record':q.name.replace('__','|'),'completed':(q/'status.txt').exists(),'error':(q/'error.txt').read_text() if (q/'error.txt').exists() else ''}
 if (q/'communications.tsv').exists():
  t=pd.read_csv(q/'communications.tsv',sep='\t');row['significant_edges']=len(t);t['source_record']=row['source_record'];tabs.append(t)
 rows.append(row)
pd.DataFrame(rows).to_csv(O/'record_execution.tsv',sep='\t',index=False)
if tabs:
 t=pd.concat(tabs);t=t.merge(audit[['source_record','F','study']],on='source_record',validate='many_to_one');t.to_csv(O/'all_record_state_specific_communications.tsv.gz',sep='\t',index=False);t.groupby(['source','target','ligand','receptor','F'],observed=True).agg(source_records=('source_record','nunique'),mean_probability=('prob','mean')).reset_index().to_csv(O/'stage_state_edge_recurrence.tsv',sep='\t',index=False)
v={'status':'completed' if rows and all(r['completed'] for r in rows) else 'partial','R_returncode':z.returncode,'eligible_records':len(rows),'completed_records':sum(r['completed'] for r in rows),'all_source_records':55,'identity':'source-scopedrecords, crossstudyuniquepersonsnotproven','method':'Original installed CellChat HumanDB, perrecord full supported cells, numeric myeloid states, min10cells/group. No inferred clinical causality. Ineligible records preserved in inputaudit.'};(O/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
