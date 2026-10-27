from pathlib import Path
import subprocess,json,sys,traceback,datetime
P=Path(__file__).parent;O=P/'F5_setup_v1';V=O/'pyscenic_venv';assert not V.exists()
pins=['numpy==1.23.5','pandas==1.5.3','scipy==1.10.1','pyarrow==14.0.2','numba==0.57.1','scikit-learn==1.3.2','dask==2023.5.0','distributed==2023.5.0','pyscenic==0.12.1','arboreto==0.1.6','ctxcore==0.2.0','loompy==3.0.8','umap-learn==0.5.5','anndata==0.9.2','matplotlib==3.7.5','networkx==3.1']
(O/'requested_versions.json').write_text(json.dumps(pins,indent=2))
try:
 with open(O/'installation.log','w') as log:
  subprocess.run([sys.executable,'-m','venv',str(V)],stdout=log,stderr=subprocess.STDOUT,check=True)
  subprocess.run([str(V/'bin/python'),'-m','pip','install','--disable-pip-version-check',*pins],stdout=log,stderr=subprocess.STDOUT,check=True)
  subprocess.run([str(V/'bin/python'),'-m','pip','check'],stdout=log,stderr=subprocess.STDOUT,check=True)
  subprocess.run([str(V/'bin/pyscenic'),'--help'],stdout=log,stderr=subprocess.STDOUT,check=True)
 r=subprocess.run([str(V/'bin/python'),'-m','pip','freeze'],capture_output=True,text=True,check=True);(O/'pip_freeze.txt').write_text(r.stdout);(O/'installation_complete.json').write_text(json.dumps({'status':'complete','python':str(V/'bin/python'),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
except Exception:(O/'installation_error.txt').write_text(traceback.format_exc());raise
