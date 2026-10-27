from pathlib import Path
import pandas as pd,json,subprocess,os
P=Path(__file__).parent;D=P.parent/'analysis_stage1_20261017_075157UTC';old=P/'GSE202379_core39_QC';new=P/'GSE202379_core39_QC_MT_NA_v2';new.mkdir(exist_ok=False)
genes=pd.read_csv(old/'genes.tsv',header=None).iloc[:,0].astype(str);matches=genes[genes.str.match(r'^(MT-|MT\.|MT_)',case=False)]
assert len(matches)==0
m=pd.read_csv(old/'cell_metadata.tsv',sep='\t',index_col=0);m['percent_mt']=float('nan');m['mt_QC_applied']=False;m.to_csv(new/'cell_metadata.tsv',sep='\t')
for f in ['counts.mtx','genes.tsv','per_donor_QC.tsv']:(new/f).symlink_to(old/f)
v=json.loads((old/'validation.json').read_text());v.update(status='partial',mt_QC_applied=False,mt_percent='NA, not zero',blocker='Supplied gene-symbol feature list contains no MT-/MT./MT_ features. UMI/gene3MAD executed; MT quality not identifiable from supplied filtered matrix. Requires unavailable upstream features.',supersedes='GSE202379_core39_QC percent_mt=0 and full3MAD claim; counts selection unchanged because zero-MAD MT never filtered cells.')
(new/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v),flush=True)
lib=P/'R_harmony_classic';lib.mkdir(exist_ok=False);env=os.environ.copy();n=D/'dependency_irgsea_repair_v3/native_dependencies/lib';env['LD_PRELOAD']=str(n/'libstdc++.so.6')+':'+str(n/'libdeflate.so.0');env['LD_LIBRARY_PATH']=str(n)
r=f".libPaths(c('{lib}','{P}/R_library','{D}/irGSEA_dependencies_v1/R_project_library_v1','{D}/trajectory_dependencies_v1/R_project_library_v1',.libPaths())); options(timeout=180);install.packages('https://cran.r-project.org/src/contrib/Archive/harmony/harmony_1.2.4.tar.gz',repos=NULL,type='source',lib='{lib}');stopifnot(as.character(packageVersion('harmony'))=='1.2.4');writeLines(capture.output(sessionInfo()),'{lib}/install_sessionInfo.txt')"
(P/'install_harmony_classic.R').write_text(r)
with open(P/'harmony_classic_install.log','w') as log:z=subprocess.run(['/usr/bin/Rscript','--vanilla',str(P/'install_harmony_classic.R')],env=env,stdout=log,stderr=subprocess.STDOUT,timeout=500)
assert z.returncode==0;print('Harmony classic1.2.4 installed',flush=True)
