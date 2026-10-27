from pathlib import Path
import requests,json,hashlib,subprocess,traceback,datetime,re
P=Path(__file__).parent;S=P/'F5_setup_v1';V=S/'pyscenic_venv';O=S/'compatibility_repair_v2';O.mkdir(exist_ok=True)
try:
 url='https://raw.githubusercontent.com/dask/dask/2024.4.1/dask/dataframe/accessor.py'
 replacement='def _bind_property(cls, pd_cls, attr, min_version=None):\n    def func(self):\n        return self._property_map(attr)\n    func.__name__ = attr\n    func.__qualname__ = f"{cls.__name__}.{attr}"\n    original_prop = getattr(pd_cls, attr)\n    if isinstance(original_prop, property):\n        method = original_prop.fget\n    elif isinstance(original_prop, functools.cached_property):\n        method = original_prop.func\n    else:\n        method = original_prop\n    try:\n        func.__wrapped__ = method\n    except Exception:\n        pass\n    setattr(cls, attr, property(derived_from(pd_cls, version=min_version)(func)))\n\n'
 (O/'official_function_excerpt_from_web.py').write_text(replacement)
 pattern=r'def _bind_property\(.*?(?=\ndef )'
 q=V/'lib/python3.11/site-packages/dask/dataframe/accessor.py';old=q.read_text();(O/'original_accessor_2023_5_0.py').write_text(old);assert len(re.findall(pattern,old,re.S))==1
 new=re.sub(pattern,lambda _:replacement,old,flags=re.S)
 if 'import functools\n' not in new:new=new.replace('import warnings','import functools\nimport warnings',1)
 compile(new,str(q),'exec');q.write_text(new)
 (O/'patch_audit.json').write_text(json.dumps({'official_source_url':url,'official_function_excerpt_sha256':hashlib.sha256(replacement.encode()).hexdigest(),'old_sha256':hashlib.sha256(old.encode()).hexdigest(),'new_sha256':hashlib.sha256(new.encode()).hexdigest(),'scope':'Backport only official _bind_property and functools import for Python3.11.9+ descriptor compatibility; no inference code changes','base_dask_version':'2023.5.0'},indent=2))
 with open(O/'verification.log','w') as log:
  for args in [[str(V/'bin/python'),'-m','pip','check'],[str(V/'bin/pyscenic'),'--help']]:subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,check=True)
 freeze=subprocess.check_output([str(V/'bin/python'),'-m','pip','freeze'],text=True);(S/'pip_freeze.txt').write_text(freeze)
 (S/'installation_error.txt').rename(O/'initial_installation_error.txt')
 (S/'installation_complete.json').write_text(json.dumps({'status':'complete','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'versions_file':'pip_freeze.txt','compatibility_patch':'compatibility_repair_v2/patch_audit.json','pip_check':True,'pyscenic_help':True},indent=2))
 I=P/'F5_GSE244832_regulon_v1'
 for n in ['pipeline_error.txt','pipeline_state.json','RSS_wait_or_run_error.txt']:
  q=I/n
  if q.exists():q.rename(O/n)
 launched={}
 for script in ['stage1_F5_pipeline_worker_v1.py','stage1_F5_RSS_wait_worker_v1.py']:
  with open(O/(script+'.log'),'ab',buffering=0) as log:p=subprocess.Popen([str(V/'bin/python'),str(P/script)],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,cwd=P);launched[script]=p.pid
 (O/'restarted_tasks.json').write_text(json.dumps(launched,indent=2));print(json.dumps(launched),flush=True)
except Exception:(O/'repair_error.txt').write_text(traceback.format_exc());raise
