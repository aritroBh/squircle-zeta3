import subprocess
PY='<project>/.venv/bin/python'
W='<project>/research/08_theory/WT5_weight5_all_k/'
r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_F_cobracket.py'], cwd=W, capture_output=True, text=True, timeout=3000)
print(r.stdout); print('STDERR:', r.stderr[-1500:])
