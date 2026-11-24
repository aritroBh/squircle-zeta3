"""flat_sum_check.py -- Theorem A check: partial sums sum_{n<=320} mu_{k,n} c^(2n) of the EXACT flat-corner
coefficients (flat_exact_N320_K3.json, Q[z2,z3]) against the independent series values M_k (values_*.txt)."""
import json
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 330
tab = json.load(open('flat_exact_N320_K3.json'))
Z = {'(0, 0, 0)': mp.mpf(1), '(1, 0, 0)': mp.zeta(2), '(0, 1, 0)': mp.zeta(3)}
def val(P):
    return mp.fsum(mp.mpf(Fr(v).numerator)/Fr(v).denominator*Z[k] for k, v in P.items())
mu = {k: [val(tab[str(k)][str(n)]) for n in range(1, 321)] for k in range(4)}
for nm, tag, a in [('5pi/6', '5pi_6', 5*mp.pi/6), ('3pi/4', '3pi_4', 3*mp.pi/4), ('2pi/3', '2pi_3', 2*mp.pi/3),
                   ('2.0', 'gen2p0', mp.mpf(2)), ('pi/2', 'pi_2', mp.pi/2)]:
    c2 = mp.cos(a/2)**2
    out = []
    for k in range(4):
        ref = mp.mpf([l for l in open(f'values_{tag}.txt') if l.startswith(f'M{k} =')][0].split('=')[1])
        s = mp.fsum(mu[k][n - 1]*c2**n for n in range(1, 321))
        out.append(mp.nstr(abs(s - ref), 3))
    print(f'{nm:>6} (c^2={mp.nstr(c2, 5)}, c^642 ~ {mp.nstr(c2**321, 3)}): |partial sum - M_k|, k=0..3:', out)
