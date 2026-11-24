# Planner: recompute every Li5 value ZQ5 stored for one D (planner_li5, own code) at 1280 digits and compare.
# Usage: planner_check_vals.py D PKLNAME   -> writes planner/d{D}_li5_planner.pkl (all values) and prints bad indices.
import sys, json, pickle, time
from fractions import Fraction as F
import mpmath as mp
sys.path.insert(0, '<project>/research/08_theory/ZQ5_zagier_weight5/planner')
from planner_li5 import li5
D, PK = int(sys.argv[1]), sys.argv[2]
HERE = f'<project>/research/08_theory/ZQ5_zagier_weight5/out/d{D}/'
J = json.load(open(HERE + 'candidates.json'))
raw = pickle.load(open(HERE + PK, 'rb'))
mp.mp.dps = 1300
s = mp.sqrt(D)
mine, bad = {}, []
t0 = time.time()
for i, c in enumerate(J['cands']):
    a, b = F(c['x'][0]), F(c['x'][1])
    x = mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator * s
    v = li5(x, 1270)
    mine[i] = mp.nstr(v, 1275)
    d = abs(v - mp.mpf(raw[f'li5:{i}']))
    if d > mp.mpf(10)**(-1150):
        bad.append(i); print(f'BAD {i:3d} x={float(x):+.6f} |planner-stored|={mp.nstr(d, 3)}', flush=True)
pickle.dump(mine, open(f'd{D}_li5_planner.pkl', 'wb'))
print(f'D={D} {PK}: {len(J["cands"])} values, {len(bad)} bad (> 1e-1150), {time.time()-t0:.0f}s')
print('BADLIST', ','.join(map(str, bad)))
