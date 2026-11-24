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
Li3_args16 = {"Li3(1/phi)": 1/phi, "Li3(-1/phi)": -1/phi, "Li3(1/phi^3)": 1/phi**3,
            "Li3(-1/phi^3)": -1/phi**3, "Li3(1/(2phi))": 1/(2*phi), "Li3(2/phi^2)": 2/phi**2,
            "Li3(1/(2phi^2))": 1/(2*phi**2)}
li3 = {nm: Li(3, x) for nm, x in Li3_args16.items()}
li2_13 = Li(2, 1/phi**3)
B3_NAMES16 = ["zeta3"] + list(Li3_args16) + ["pi2*L", "pi2*lam", "Li2(1/phi^3)*L", "Li2(1/phi^3)*lam",
            "L^3", "L^2*lam", "L*lam^2", "lam^3"]
B3_16 = [z3] + list(li3.values()) + [pi2*L, pi2*lam, li2_13*L, li2_13*lam, L**3, L**2*lam, L*lam**2, lam**3]
nb16 = len(B3_16)
W3 = mp.mpf(open(WD+"gol3_W3_hpl_1150.txt").read())
vals = [W3] + B3_16 + [r5*x for x in B3_16]
rows = L7.find_relation_rows(vals, 1000)
mp.mp.dps = 1120
cands = []
for r in rows:
    if r[0] == 0: continue
    g = L7.normalize(r)
    resid = sum(c*v for c, v in zip(g, vals))
    if abs(resid) < mp.mpf('1e-1000'):
        cands.append(g)
print("valid candidate relations for W3(gol3):", len(cands))
for g in cands:
    print("  a0 =", g[0], " height", max(abs(c) for c in g))
    terms = []
    for i in range(nb16):
        p, q = g[i+1], g[i+1+nb16]
        if p or q:
            terms.append(f"({Fr(-p, g[0])}{'+' if q>=0 else ''}{Fr(-q, g[0])}s5)*{B3_NAMES16[i]}")
    print("   ", " + ".join(terms))
