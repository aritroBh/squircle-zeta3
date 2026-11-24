import time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from cor9_lib import L7
from fractions import Fraction as Fr
mp.mp.dps = 1120
r5 = mp.sqrt(5); phi = (1+r5)/2
L = mp.log(phi); lam = mp.log(2); pi2 = mp.pi**2
Li = lambda s, x: mp.polylog(s, x)
z3 = mp.zeta(3)
Li3_args = {"Li3(1/phi)": 1/phi, "Li3(-1/phi)": -1/phi, "Li3(1/phi^2)": 1/phi**2,
            "Li3(1/phi^3)": 1/phi**3, "Li3(1/(2phi))": 1/(2*phi)}
li3 = {nm: Li(3, x) for nm, x in Li3_args.items()}
li2_13 = Li(2, 1/phi**3)
B3_NAMES = ["zeta3"] + list(Li3_args) + ["pi2*L", "pi2*lam", "Li2(1/phi^3)*L", "Li2(1/phi^3)*lam",
            "L^3", "L^2*lam", "L*lam^2", "lam^3"]
B3 = [z3] + list(li3.values()) + [pi2*L, pi2*lam, li2_13*L, li2_13*lam, L**3, L**2*lam, L*lam**2, lam**3]
nb = len(B3)
vals = B3 + [r5*x for x in B3]
rows = L7.find_relation_rows(vals, 900)
mp.mp.dps = 1120
print("short rows among basis alone (resid@1120 < 1e-1000):")
nrel = 0
for r in rows:
    g = L7.normalize(r)
    resid = sum(c*v for c, v in zip(g, vals))
    if abs(resid) < mp.mpf('1e-1000'):
        nrel += 1
        terms = []
        for i in range(nb):
            p, q = g[i], g[i+nb]
            if p or q:
                terms.append(f"({p}{'+' if q>=0 else ''}{q}*s5)*{B3_NAMES[i]}")
        print(f"  height {max(abs(c) for c in g)}: {' + '.join(terms)}")
print("total pure-basis relations:", nrel)
