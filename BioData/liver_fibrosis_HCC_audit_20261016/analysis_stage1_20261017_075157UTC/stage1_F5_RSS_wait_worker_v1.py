from pathlib import Path
import time,subprocess,json,traceback,datetime
P=Path(__file__).parent;I=P/'F5_GSE244832_regulon_v1'
try:
 start=time.monotonic()
 while not (I/'CSI_review_v1/validation.json').exists():
  if (I/'pipeline_error.txt').exists():raise RuntimeError('Main pipeline failed; RSS not started. See pipeline_error.txt')
  if time.monotonic()-start>21600:raise TimeoutError('Waiting six hours; RSS not started')
  time.sleep(30)
 with open(I/'RSS_postprocess.log','w') as f:r=subprocess.run([str(P/'F5_setup_v1/pyscenic_venv/bin/python'),str(P/'stage1_F5_RSS_postprocess_v1.py')],stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0
except Exception:(I/'RSS_wait_or_run_error.txt').write_text(traceback.format_exc());raise
