import time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
from cor9_lib import L7
from fractions import Fraction as Fr
mp.mp.dps = 1120
W3 = mp.mpf(open(WD+"b2sq_W3_hpl_1150.txt").read())
b = mp.sqrt(2)-1; rt2 = mp.sqrt(2)
Llog = mp.log(1+rt2); lam = mp.log(2)
Li = lambda s, x: mp.polylog(s, x)
z3, pi2 = mp.zeta(3), mp.pi**2
B3_NAMES = ["zeta3", "Li3(b)", "Li3(-b)", "Li3(1/sqrt2)", "L*Li2(b)", "lam*Li2(b)",
            "pi2*L", "pi2*lam", "lam^3", "lam^2*L", "lam*L^2", "L^3"]
Li2b = Li(2, b)
B3 = [z3, Li(3, b), Li(3, -b), Li(3, 1/rt2), Llog*Li2b, lam*Li2b,
      pi2*Llog, pi2*lam, lam**3, lam**2*Llog, lam*Llog**2, Llog**3]
vals = [W3] + B3 + [rt2*x for x in B3]
cands = {}
for dps in (700, 1000):
    t0=time.time()
    rows = L7.find_relation_rows(vals, dps)
    mp.mp.dps = 1120
    # pick the first row with nonzero target coefficient
    row = next(r for r in rows if r[0] != 0)
    cands[dps] = L7.normalize(row)
    resid = sum(c*v for c, v in zip(cands[dps], vals))
    print(f"dps {dps}: LLL {time.time()-t0:.0f}s, height {max(abs(c) for c in cands[dps])}, resid@1120 {mp.nstr(abs(resid),3)}", flush=True)
print("identical:", cands[700] == cands[1000])
a = cands[1000]; a0 = a[0]
print(f"\nW3((sqrt2-1)^2), a0 = {a0};  W3 = sum (coeff in Q(sqrt2)) * basis:")
for i in range(1, 13):
    p, q = -a[i], -a[i+12]
    if p or q:
        print(f"  ({Fr(p, a0)} {'+' if q>=0 else ''}{Fr(q, a0)}*sqrt2) * {B3_NAMES[i-1]}")
json.dump({"a0": a0, "coeffs": a[1:]}, open(WD+"b2sq_W3_relation.json", "w"))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_27: b2sq W3 LLL dps 700/1000 identical={cands[700]==cands[1000]}.\n")
