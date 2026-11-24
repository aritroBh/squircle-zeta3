# Planner: recompute ZQ5's stored D=21 Li5 values (vals_b10.pkl) by direct series at 1250 digits,
# starting with the support of kernel combo K15 (flagged ~ -9.84e-726 by ZQ5). Prints per candidate.
import json, pickle
from fractions import Fraction as F
import mpmath as mp
HERE = '<project>/research/08_theory/ZQ5_zagier_weight5/out/d21/'
J = json.load(open(HERE + 'candidates.json'))
K = json.load(open(HERE + 'kernel_b10.json'))['K']
raw = pickle.load(open(HERE + 'vals_b10.pkl', 'rb'))
mp.mp.dps = 1280
s = mp.sqrt(J['D'])
def xval(c):
    a, b = F(c['x'][0]), F(c['x'][1])
    return mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator * s
def li5(x):
    tot = mp.mpf(0); p = x; n = 1; eps = mp.mpf(10)**(-1265)
    while abs(p) > eps:
        tot += p / mp.mpf(n)**5; n += 1; p *= x
    return tot
supp = [i for i, c in enumerate(K[15]) if c]
order = supp + [i for i in range(len(J['cands'])) if i not in supp]
print('K15 support:', supp, flush=True)
mine = {}
for i in order:
    x = xval(J['cands'][i])
    v = li5(x); mine[i] = v
    d = abs(v - mp.mpf(raw[f'li5:{i}']))
    print(f'{i:3d} x={float(x):+.6f} |mine-stored|={mp.nstr(d, 3)}', flush=True)
    if i == supp[-1]:
        k15 = sum(K[15][j] * mine[j] for j in supp)
        print('K15 with my values =', mp.nstr(k15, 5), flush=True)
