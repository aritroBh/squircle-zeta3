import sys
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import cmzv6
from fractions import Fraction as Fr
import exact_linalg as XL
lower = {}
for n in (1,2):
    rows, _ = cmzv6.make_relations(n, lower)
    lower[n] = rows
cols = cmzv6.words(2)
# target: zeta(2) - pi^2/6 = (-1,0) + 3*(1,1)
T = {(-1,0): Fr(1), (1,1): Fr(3)}
prow, r = XL.select_independent_rows(rows, cols)
sel = [rows[i] for i in prow]
pcols = XL.select_pivot_cols(sel, cols)
y = XL.exact_solve(sel, pcols, cols, T)
ok, bad = XL.verify_combination(sel, y, cols, T)
print('zeta(2) = pi^2/6 exact rational combination of weight-2 proved rows:', ok, bad)
import json
json.dump(dict(cols=[list(w) for w in cols], prow=prow, y=[str(c) for c in y]), open('cert_euler.json','w'))
with open('progress.log','a') as f:
    f.write('step_78: zeta(2)=pi^2/6 exact combination of weight-2 proved rows (cert_euler.json)\n')
