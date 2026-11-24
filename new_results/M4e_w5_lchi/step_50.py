import subprocess, time, json
M4E = '<project>/research/08_theory/M4e_w5_lchi'
PY = '<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice', PY, 'lchi_final.py', '1500'], cwd=M4E, capture_output=True, text=True, timeout=2700)
print(r.stdout); print("STDERR:", r.stderr[-1500:]); print(f"elapsed {time.monotonic()-t0:.0f}s")
# compare the two relations
from fractions import Fraction as F
d12 = json.load(open(M4E+'/lchi_final_1200.json')); d15 = json.load(open(M4E+'/lchi_final_1500.json'))
print("c0:", d12['c0'], d15['c0'])
same = (d12['ca']==d15['ca'] and d12['cb']==d15['cb'] and d12['aL']==d15['aL'] and d12['bL']==d15['bL'])
print("atom coefficients identical at 1200 and 1500:", same)
open(M4E+'/progress.log','a').write("step_47: lchi_final.py 1500 (stability check)\n")
