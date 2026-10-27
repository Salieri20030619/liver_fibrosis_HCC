from pathlib import Path
import time,json,subprocess
p=Path(__file__).parent
for i in range(240):
 s=json.loads((p/'F4_irGSEA_official_v2/state.json').read_text())
 if s['status']=='R_steps_completed':
  raise SystemExit(subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python',str(p/'stage1_irGSEA_validate_plot_v3.py')]).returncode)
 if s['status']=='blocked':raise RuntimeError('Original method task blocked; no plots created')
 time.sleep(15)
raise TimeoutError('No completed method results after 1h')
