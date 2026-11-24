import shutil, os, subprocess, time
base = '<project>'
work = f'{base}/research/08_theory/M5a_next_coefficient'
os.chdir(work)
# Copy the needed scripts from paper (read-only there; copy into my folder)
for f in ['ibp_reduce.py', 'mk_integrands.py', 'qsqrt2.py']:
    shutil.copy(f'{base}/paper/squircle_note/scripts/{f}', f'{work}/{f}')
print("copied:", os.listdir(work))
PY = f'{base}/.venv/bin/python'
# First quick check: time the integrand generation for k=5
t0 = time.time()
r = subprocess.run([PY, '-B', '-c', '''
import time, sys
sys.dont_write_bytecode=True
t0=time.time()
from mk_integrands import v_integrand
R5 = v_integrand(5)
print("v_integrand(5) built in", time.time()-t0, "s; keys:", sorted(R5.keys()))
'''], capture_output=True, text=True, cwd=work)
print(r.stdout, r.stderr)
print("elapsed:", time.time()-t0)