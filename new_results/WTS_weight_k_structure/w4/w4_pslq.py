"""EXPLORATORY identification (not a proof): try to express c_2, delta, c_1, c_0 by PSLQ over small weight-<=4 bases.
Inputs: results/poly_w4_dps120_k4-10.json (fit through k=4,5,6) and results/w4_k3_dps120.json.  dps 120."""
import json
from w4_common import *
set_dps(120); print("mp.mp.dps =", mp.mp.dps)
P = json.load(open(os.path.join(RES, "poly_w4_dps120_k4-10.json")))
c = [mp.mpf(x) for x in P["coeffs"]]
W43 = mp.mpf(json.load(open(os.path.join(RES, "w4_k3_dps120.json")))["value"])
delta = c[0] + 3*c[1] + 9*c[2] - W43     # P(3) - W_4(3)
print("delta = P(3) - W_4(3) =", mp.nstr(delta, 50))
lam = mp.log(2); z3 = mp.zeta(3); pi2 = mp.pi**2; L = mp.log(1 + mp.sqrt(2)); r2 = mp.sqrt(2)
bases = {
 "Q4 (pi^4, pi^2 lam^2, lam^4, z3 lam, Li4(1/2))": ["pi^4", "pi^2 lam^2", "lam^4", "z3 lam", "Li4(1/2)"],
}
vals = {"pi^4": pi2**2, "pi^2 lam^2": pi2*lam**2, "lam^4": lam**4, "z3 lam": z3*lam, "Li4(1/2)": mp.polylog(4, mp.mpf(1)/2)}
for tname, target in (("c_2", c[2]), ("delta", delta), ("c_1", c[1]), ("c_0", c[0])):
    for bname, names in bases.items():
        vec = [target] + [vals[n] for n in names]
        rel = mp.pslq(vec, maxcoeff=10**6, maxsteps=10**6, tol=mp.mpf(10)**(-100))
        print(f"{tname} over {bname}: relation =", rel)
        if rel:
            res = sum(r*v for r, v in zip(rel, vec))
            print("   residual", mp.nstr(res, 5), " (dps", mp.mp.dps, ")")
