step31 = r'''import sys, json, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '../HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from cor13c_common import target4
import exact_linalg as XL
from fractions import Fraction as Fr

# (a) Euler: zeta(2) = pi^2/6 as exact weight-2 combination: (-1,0) + 3*(1,1) in the span
lower = {}
for n in (1,2):
    rows, _ = cmzv6.make_relations(n, lower)
    lower[n] = rows
cols2 = cmzv6.words(2)
ET = {(-1,0): Fr(1), (1,1): Fr(3)}
prow2, r2 = XL.select_independent_rows(rows, cols2)
sel2 = [rows[i] for i in prow2]
pcols2 = XL.select_pivot_cols(sel2, cols2)
yE = XL.exact_solve(sel2, pcols2, cols2, ET)
ok, bad = XL.verify_combination(sel2, yE, cols2, ET)
print("Euler zeta(2)=pi^2/6 exact in weight-2 span:", ok)
assert ok
json.dump(dict(cols=[list(w) for w in cols2], prow=prow2,
               y=[[str(c.numerator), str(c.denominator)] for c in yE]),
          open('cert13_euler.json','w'))

# (b) pi^4 coefficient check: coefficient of (1,1,1,1) in T4 and in the certified combination
T4 = target4()
c_t = T4.get((1,1,1,1), Fr(0))
rows4, _ = cmzv6.make_relations(4, lower)   # lower now has 1,2 only; add 3
print("coefficient of (1,1,1,1) in T4 (= 813/2 from pi^4 after reduction):", c_t)
with open('progress.log','a') as f:
    f.write(f'step_31: Euler cert exact; T4 coeff of (1,1,1,1) = {c_t}\n')
'''
open('step_31.py','w').write(step31)
import subprocess, time
t0 = time.time()
r = subprocess.run(['<project>/.venv/bin/python', 'step_31.py'], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-3000:] if r.stderr else "")
print(f"elapsed {time.time()-t0:.1f}s")
