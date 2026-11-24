import subprocess
PY='<project>/.venv/bin/python'
W='<project>/research/08_theory/WT5_weight5_all_k/'
for t in ['e1','e2','e3']:
    r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_precision.py', t, '600'], cwd=W, capture_output=True, text=True, timeout=1500)
    print(r.stdout[-1200:]); print('STDERR:', r.stderr[-500:])
