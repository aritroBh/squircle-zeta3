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
Li3_args = {"Li3(1/phi)": 1/phi, "Li3(1/phi^3)": 1/phi**3, "Li3(1/(2phi))": 1/(2*phi), "Li3(2/phi^2)": 2/phi**2}
li3 = {nm: Li(3, x) for nm, x in Li3_args.items()}
li2_13 = Li(2, 1/phi**3)
B3_NAMES = ["zeta3"] + list(Li3_args) + ["pi2*L", "pi2*lam", "Li2(1/phi^3)*L", "Li2(1/phi^3)*lam",
            "L^3", "L^2*lam", "L*lam^2", "lam^3"]
B3 = [z3] + list(li3.values()) + [pi2*L, pi2*lam, li2_13*L, li2_13*lam, L**3, L**2*lam, L*lam**2, lam**3]
nb = len(B3)
vals = B3 + [r5*x for x in B3]
rows = L7.find_relation_rows(vals, 900)
mp.mp.dps = 1120
pure = sum(1 for r in rows if abs(sum(c*v for c, v in zip(L7.normalize(r), vals))) < mp.mpf('1e-1000'))
print("pure-basis relations in 13-el basis:", pure, "(basis size", nb, ")")
for tag in ("gol1", "gol3"):
    W3 = mp.mpf(open(WD+f"{tag}_W3_hpl_1150.txt").read())
    vals2 = [W3] + B3 + [r5*x for x in B3]
    cands = {}
    for dps in (700, 1000):
        rows = L7.find_relation_rows(vals2, dps)
        mp.mp.dps = 1120
        row = next(r for r in rows if r[0] != 0)
        cands[dps] = L7.normalize(row)
        resid = sum(c*v for c, v in zip(cands[dps], vals2))
        print(f"{tag} dps {dps}: height {max(abs(c) for c in cands[dps])}, resid@1120 {mp.nstr(abs(resid),3)}", flush=True)
    print(tag, "identical:", cands[700] == cands[1000])
    if cands[700] == cands[1000]:
        a = cands[1000]; a0 = a[0]
        print(f"W3({tag}), a0 = {a0}:")
        for i in range(1, nb+1):
            p, q = -a[i], -a[i+nb]
            if p or q:
                print(f"  ({Fr(p, a0)} {'+' if q>=0 else ''}{Fr(q, a0)}*sqrt5) * {B3_NAMES[i-1]}")
        json.dump({"a0": a0, "coeffs": a[1:], "basis": B3_NAMES}, open(WD+f"{tag}_W3_relation.json", "w"))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_41: golden W3 LLL with independent 13-el basis.\n")
