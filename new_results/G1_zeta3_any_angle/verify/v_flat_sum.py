"""v_flat_sum.py -- Theorem A series sum_{n<=320} mu_{k,n} c^(2n) (solver's exact table, z2,z3 substituted)
vs the verifier's independent Cauchy-contour M_k (runB, ~1e-46), at angles incl. fresh ones."""
import json, mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 60
tab = json.load(open('flat_exact_N320_K3.json'))
z2, z3 = mp.zeta(2), mp.zeta(3)
mu = {k: [mp.mpf(0)] + [sum(mp.mpf(Fr(v).numerator)/Fr(v).denominator*z2**eval(key)[0]*z3**eval(key)[1]
                         for key, v in tab[str(k)][str(n)].items()) for n in range(1, 321)] for k in range(4)}
def rd(fn):
    return {int(l.split('=')[0].strip()[1]): mp.mpf(l.split('=')[1]) for l in open(fn) if l.startswith('M')}
for tag, a in [('new3p0', mp.mpf(3)), ('new2p6', mp.mpf(26)/10), ('5pi_6', 5*mp.pi/6), ('3pi_4', 3*mp.pi/4),
               ('2pi_3', 2*mp.pi/3), ('gen2p0', mp.mpf(2)), ('pi_2', mp.pi/2)]:
    c2 = mp.cos(a/2)**2; V = rd(f'runB/vc_{tag}.txt')
    print(f'{tag:>7} c^2={mp.nstr(c2,5)} c^642={mp.nstr(c2**321,3)}: ' +
          ' '.join(f'|sum-M{k}|={mp.nstr(abs(mp.fsum(mu[k][n]*c2**n for n in range(1, 321)) - V[k]), 3)}' for k in range(4)))
