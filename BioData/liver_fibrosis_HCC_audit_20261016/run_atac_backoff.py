import pathlib,subprocess,json,datetime
r=pathlib.Path('/mnt/BioData/bio/wangmingyuan/liver_fibrosis_HCC_audit_20261016')
jobs=[('GSE227265','GSE227265_BinarizedPeaks.mtx.gz',369682663),('GSE227265_TF','GSE227265_TFMotifActivity.csv.gz',70534594)]
for g,n,size in jobs:
 if (r/(n+'.complete')).exists():continue
 url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE227nnn/GSE227265/suppl/'+n
 with (r/(g+'_low_concurrency_'+datetime.datetime.now().strftime('%H%M%S')+'.log')).open('w') as log:
  subprocess.run(['python3',str(r/'resume_backoff.py'),g,n,str(size),str(r),str(r/(g+'.absent_prefix')),url],stdout=log,stderr=subprocess.STDOUT)
 if (r/(n+'.complete')).exists() and g=='GSE227265':
  subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python',str(r/'validate_peak.py')],check=True)
print('Bounded two-connection ATAC retry round finished',flush=True)
