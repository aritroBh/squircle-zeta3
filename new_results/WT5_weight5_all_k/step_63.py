import subprocess
PY='<project>/.venv/bin/python'
W='<project>/research/08_theory/WT5_weight5_all_k/'
for i in ['0','1']:
    r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_I_assemble.py', 'e3', i], cwd=W, capture_output=True, text=True, timeout=600)
    print('===== e3 rel', i); print(r.stdout); print('STDERR:', r.stderr[-300:])
