import subprocess
PY='<project>/.venv/bin/python'
W='<project>/research/08_theory/WT5_weight5_all_k/'
r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_precision.py', 'e2', '5400', '5800'], cwd=W, capture_output=True, text=True, timeout=3000)
print(r.stdout[-900:]); print('STDERR:', r.stderr[-300:])
