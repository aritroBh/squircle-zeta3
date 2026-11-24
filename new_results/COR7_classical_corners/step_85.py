import sys, json
sys.path.insert(0, '.')
import mpmath as mp
import cor7_lib as L
from fractions import Fraction as Fr

mp.mp.dps = 1120
W3 = mp.mpf(open("b21_W3_hpl_1150.txt").read())
b = mp.sqrt(2)-1
rt2 = mp.sqrt(2)
Llog = mp.log(1+mp.sqrt(2))   # L = log(1+sqrt2)
lam = mp.log(2)
def Li(s, x): return mp.polylog(s, x)
z3, pi2 = mp.zeta(3), mp.pi**2
B3_NAMES = ["zeta3", "Li3(b)", "Li3(-b)", "Li3(1/sqrt2)", "L*Li2(b)", "lam*Li2(b)",
            "pi2*L", "pi2*lam", "lam^3", "lam^2*L", "lam*L^2", "L^3"]
Li2b = Li(2, b)
B3 = [z3, Li(3, b), Li(3, -b), Li(3, 1/rt2), Llog*Li2b, lam*Li2b,
      pi2*Llog, pi2*lam, lam**3, lam**2*Llog, lam*Llog**2, Llog**3]
vals = [W3] + B3 + [rt2*x for x in B3]
cands = {}
for dps in (700, 1000):
    rows = L.find_relation_rows(vals, dps)
    cands[dps] = L.normalize(rows[0])
    mp.mp.dps = 1120
    resid = sum(c*v for c, v in zip(cands[dps], vals))
    print(f"dps {dps}: height {max(abs(c) for c in cands[dps])}, resid@1120 {mp.nstr(abs(resid),3)}", flush=True)
print("identical:", cands[700] == cands[1000])
a = cands[1000]; a0 = a[0]
print("\nW3(sqrt2-1) =  (coeffs in Q(sqrt2): c = p + q*sqrt2)")
for i in range(1, 13):
    p, q = -a[i], -a[i+12]
    if p or q:
        print(f"  ({Fr(p, a0)} {'+' if q>=0 else ''}{Fr(q, a0)}*sqrt2) * {B3_NAMES[i-1]}")
json.dump({"a0": a0, "coeffs": a[1:]}, open("b21_W3_relation.json", "w"))
