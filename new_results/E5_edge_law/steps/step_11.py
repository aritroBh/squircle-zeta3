# step_11.py: PSLQ of the independent numeric P^(5)(1)/sqrt5 (rhombus 1,2, step_10) on
# {1,z2,z3,z2^2,z2z3,z5}; the zeta(5) coordinate must equal sum_e Z_e^(5)/sqrt5 = 1152/25.
import mpmath as mp, json
from fractions import Fraction
mp.mp.dps = 100
P5 = mp.mpf(json.load(open('steps/step_10_num.json'))['P5'])
z2, z3, z5 = mp.zeta(2), mp.zeta(3), mp.zeta(5)
V = P5/mp.sqrt(5)
basis = [mp.mpf(1), z2, z3, z2**2, z2*z3, z5]
bn = ['1','z2','z3','z2^2','z2z3','z5']
rel = None
for tol in ['1e-55', '1e-50', '1e-45']:
    rel = mp.pslq(mp.matrix([V]+basis), tol=mp.mpf(tol), maxcoeff=10**8, maxsteps=10000)
    if rel is not None: break
print("relation:", rel)
q = [Fraction(int(rel[i+1]), int(-rel[0])) for i in range(6)]
print("P^(5)(1)/sqrt5 =", {bn[i]: str(q[i]) for i in range(6) if q[i] != 0})
resid = abs(V + sum(mp.mpf(rel[i+1])*basis[i]/rel[0] for i in range(6)))
print("residual:", mp.nstr(resid, 3))
exact_vec = [0, sp.Rational(0)] if False else None
import sympy as sp
exact = {'z2': sp.Rational(-1536,5), 'z3': sp.Rational(2304,5), 'z2^2': sp.Rational(-4608,125),
         'z2z3': sp.Rational(-768,25), 'z5': sp.Rational(1152,25), '1': sp.Rational(0)}
got = {bn[i]: sp.Rational(q[i]) for i in range(6)}
print("matches exact vector (incl. z5 = 1152/25, z2z3 = -768/25):",
      all(got[k] == exact[k] for k in exact))
print("zeta(5) coordinate = sqrt5 * %s = %s = sum_e Z_e^(5)" % (q[5], sp.sqrt(5)*sp.Rational(q[5])))
