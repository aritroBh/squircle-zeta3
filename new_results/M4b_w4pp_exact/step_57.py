import os, subprocess, time
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
PY = os.path.join(ROOT, '.venv/bin/python')
if os.path.exists('verifier.log'): os.remove('verifier.log')
t = time.time()
r = subprocess.run([PY, 'verify_W4pp.py'], capture_output=True, text=True, cwd=WORK)
print('wall', round(time.time()-t,1),'s  EXIT', r.returncode)
print(r.stdout)
print(r.stderr[-3000:] if r.returncode else '')
