import time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from cor9_lib import L7
mp.mp.dps = 300
r5 = mp.sqrt(5); phi = (1+r5)/2
L = mp.log(phi); lam = mp.log(2); pi2 = mp.pi**2
Li = lambda s, x: mp.polylog(s, x)
v = Li(2, 1/phi**3)
rows = L7.find_relation_rows([v, pi2, L**2, L*lam, lam**2], 250)
r = next(rr for rr in rows if rr[0] != 0)
from fractions import Fraction as Fr
g = L7.normalize(r)
resid = sum(c*x for c, x in zip(g, [v, pi2, L**2, L*lam, lam**2]))
print("Li2(1/phi^3) relation coeffs:", g, " resid:", mp.nstr(abs(resid), 3))
print("check height:", max(abs(c) for c in g))
