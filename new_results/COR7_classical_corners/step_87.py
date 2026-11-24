import sys, json, time
sys.path.insert(0, '.')
import mpmath as mp
import cor7_lib as L
from fractions import Fraction as Fr

mp.mp.dps = 1120
W4 = mp.mpf(open("b21_W4_hpl_1150.txt").read())
b = mp.sqrt(2)-1
rt2 = mp.sqrt(2)
Llog = mp.log(1+rt2)
lam = mp.log(2)
def Li(s, x): return mp.polylog(s, x)
z3, pi2, pi4 = mp.zeta(3), mp.pi**2, mp.pi**4
Li3b, Li3mb, Li3s = Li(3, b), Li(3, -b), Li(3, 1/rt2)
Li2b = Li(2, b)
B4_NAMES = ["Li4(b)", "Li4(-b)", "Li4(-b^2)", "Li4(1/sqrt2)", "Li4(1/2)", "Li4(1-1/sqrt2)", "Li4(2-sqrt2)",
            "L*Li3(b)", "L*Li3(-b)", "L*Li3(1/sqrt2)", "L*zeta3",
            "lam*Li3(b)", "lam*Li3(-b)", "lam*Li3(1/sqrt2)", "lam*zeta3",
            "Li2(b)^2", "Li2(b)*L^2", "Li2(b)*lam*L", "Li2(b)*lam^2", "pi2*Li2(b)",
            "pi4", "pi2*L^2", "pi2*lam*L", "pi2*lam^2",
            "L^4", "L^3*lam", "L^2*lam^2", "L*lam^3", "lam^4"]
B4 = [Li(4, b), Li(4, -b), Li(4, -b*b), Li(4, 1/rt2), Li(4, mp.mpf(1)/2), Li(4, 1-1/rt2), Li(4, 2-rt2),
      Llog*Li3b, Llog*Li3mb, Llog*Li3s, Llog*z3,
      lam*Li3b, lam*Li3mb, lam*Li3s, lam*z3,
      Li2b**2, Li2b*Llog**2, Li2b*lam*Llog, Li2b*lam**2, pi2*Li2b,
      pi4, pi2*Llog**2, pi2*lam*Llog, pi2*lam**2,
      Llog**4, Llog**3*lam, Llog**2*lam**2, Llog*lam**3, lam**4]
print("basis size:", len(B4))
vals = [W4] + B4 + [rt2*x for x in B4]
print("lattice dim:", len(vals))
cands = {}
for dps in (700, 1000):
    t0 = time.time()
    rows = L.find_relation_rows(vals, dps)
    cands[dps] = L.normalize(rows[0])
    mp.mp.dps = 1120
    resid = sum(c*v for c, v in zip(cands[dps], vals))
    print(f"dps {dps}: LLL {time.time()-t0:.0f}s, height {max(abs(c) for c in cands[dps])}, resid@1120 {mp.nstr(abs(resid),3)}", flush=True)
print("identical:", cands[700] == cands[1000])
a = cands[1000]; a0 = a[0]
nb = len(B4)
print("\nW4(sqrt2-1) =")
for i in range(1, nb+1):
    p, q = -a[i], -a[i+nb]
    if p or q:
        print(f"  ({Fr(p, a0)} {'+' if q>=0 else ''}{Fr(q, a0)}*sqrt2) * {B4_NAMES[i-1]}")
json.dump({"a0": a0, "coeffs": a[1:], "basis": B4_NAMES}, open("b21_W4_relation.json", "w"))
