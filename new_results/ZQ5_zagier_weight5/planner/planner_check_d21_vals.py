# Planner: independent recomputation of ZQ5's D=21 Li5 values (vals_b10.pkl) at 1250 digits.
import json, pickle, sys
from fractions import Fraction as F
import mpmath as mp
HERE = '<project>/research/08_theory/ZQ5_zagier_weight5/out/d21/'
J = json.load(open(HERE + 'candidates.json'))
raw = pickle.load(open(HERE + 'vals_b10.pkl', 'rb'))
mp.mp.dps = 1300
D = J['D']; s = mp.sqrt(D)
def li5_series(x):
    # direct series, enough terms for 1250 digits
    tot = mp.mpf(0); p = x; n = 1; eps = mp.mpf(10)**(-1260)
    while True:
        t = p / mp.mpf(n)**5
        tot += t
        if abs(p) < eps: break
        n += 1; p *= x
    return tot
bad = []
worst = mp.mpf(0)
for i, c in enumerate(J['cands']):
    a, b = F(c['x'][0]), F(c['x'][1])
    x = mp.mpf(a.numerator)/a.denominator + mp.mpf(b.numerator)/b.denominator * s
    v_me = li5_series(x)
    v_pp = mp.polylog(5, x)
    v_k = mp.mpf(raw[f'li5:{i}'])
    d1 = abs(v_me - v_k); d2 = abs(v_me - v_pp)
    worst = max(worst, d1)
    if d1 > mp.mpf(10)**(-1140):
        bad.append((i, c['x'], float(x), mp.nstr(d1, 3), mp.nstr(d2, 3)))
print('D', D, 'candidates', len(J['cands']), 'stored keys li5', sum(k.startswith('li5:') for k in raw))
print('max |mine - stored| =', mp.nstr(worst, 3))
print('values off by more than 1e-1140:', len(bad))
for t in bad: print('  ', t)
# other stored atoms
for k in ('z5', 'L', 'T', 'pi2'):
    if k in raw: print(k, 'stored digits ok?')
