from pathlib import Path
import subprocess,sys,os,json
P=Path(__file__).parent;O=P/'python_deps';O.mkdir(exist_ok=False)
cmd=[sys.executable,'-m','pip','install','--target',str(O),'--no-deps','--disable-pip-version-check','xgboost-cpu==2.1.4','shap==0.46.0','slicer==0.0.8']
r=subprocess.run(cmd,capture_output=True,text=True,timeout=180)
(P/'python_dependencies_install.log').write_text(r.stdout+'\n'+r.stderr)
print(json.dumps({'command':cmd,'returncode':r.returncode,'tail':(r.stdout+r.stderr)[-2500:]}));sys.exit(r.returncode)
