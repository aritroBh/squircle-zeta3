import subprocess
PY='<project>/.venv/bin/python'
W='<project>/research/08_theory/WT5_weight5_all_k/'
r = subprocess.run(['nice','-n','10', PY, '-B', 'wt5_K_verify.py', 'e1'], cwd=W, capture_output=True, text=True, timeout=3000)
# filter out the noisy relation prints
out=[l for l in r.stdout.splitlines() if 'relation multiplier' not in l and 'exact zero relation' not in l]
print('\n'.join(out)); print('STDERR:', r.stderr[-500:])
