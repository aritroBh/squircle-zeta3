"""w6 identification step 6 (structure guide): PSLQ of the lower-weight polynomial coefficients (lower_coef.json, dps 150 data)
against pure level-2 bases: weight 4 {z3 lam, pi^4, pi^2 lam^2, lam^4, Li4(1/2)}; weight 5 8-element basis.  Exploratory."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
lc = json.load(open(os.path.join(W, "lower_coef.json")))
mp.mp.dps = 140
print("mp.mp.dps =", mp.mp.dps)
lam = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi ** 2; h = mp.mpf(1) / 2
Li4 = mp.polylog(4, h); Li5 = mp.polylog(5, h)
B4 = [("z3 lam", z3 * lam), ("pi^4", pi2 ** 2), ("pi^2 lam^2", pi2 * lam ** 2), ("lam^4", lam ** 4), ("Li4(1/2)", Li4)]
B5 = [("z5", z5), ("pi^2 z3", pi2 * z3), ("z3 lam^2", z3 * lam ** 2), ("pi^4 lam", pi2 ** 2 * lam), ("pi^2 lam^3", pi2 * lam ** 3),
      ("lam^5", lam ** 5), ("lam Li4(1/2)", lam * Li4), ("Li5(1/2)", Li5)]
D = 120
for wt, B in ((4, B4), (5, B5)):
    cs = lc[f"W{wt}"]
    for j in range(len(cs) - 1, -1, -1):
        mp.mp.dps = 140
        c = mp.mpf(cs[j])
        rel, res, tol = pslq_try(c, B, D, maxcoeff=10 ** 9)
        print(f"weight {wt}, c_{j} = {mp.nstr(c, 25)}")
        show(f"pure level-2 basis ({len(B)})", B, rel, res, tol, D)
