import time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
mp.mp.dps = 300
r5 = mp.sqrt(5); phi = (1+r5)/2
L = mp.log(phi); lam = mp.log(2); pi2 = mp.pi**2
Li = lambda s, x: mp.polylog(s, x)
# elementary Li2 checks
for nm, x, cf in [("Li2(1/phi)", 1/phi, pi2/10 - L**2),
                  ("Li2(1/phi^2)", 1/phi**2, pi2/15 - L**2),
                  ("Li2(-1/phi)", -1/phi, -pi2/15 + L**2/2),
                  ("Li2(-1/phi^2)", -1/phi**2, -pi2/10 + L**2),
                  ("Li2(1/phi^3)", 1/phi**3, None)]:
    v = Li(2, x)
    if cf is not None:
        print(nm, "|Li2 - elementary| =", mp.nstr(abs(v-cf), 3))
    else:
        print(nm, "=", mp.nstr(v, 20))
