import subprocess
PY='<project>/.venv/bin/python'
W='<project>/research/08_theory/WT5_weight5_all_k/'
r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_H_lll.py', 'e1', '1300'], cwd=W, capture_output=True, text=True, timeout=3000)
print(r.stdout[-4000:]); print('STDERR:', r.stderr[-500:])
