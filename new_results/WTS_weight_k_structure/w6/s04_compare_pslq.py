"""w6 step 4: (a) 70-dps vs 110-dps agreement of the fitted coefficients (reliable digits);
(b) exploratory PSLQ of the leading coefficient c_4 (and the k=5 defect) against weight-6 constants.
All printed numbers come from this script (110-dps run is the high-precision source).
"""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp

a70 = json.load(open(os.path.join(W, "analysis_vals.json")))
a110 = json.load(open(os.path.join(W, "analysis_vals_110.json")))
mp.mp.dps = 110       # set after imports
print("mp.mp.dps =", mp.mp.dps)
c70 = [mp.mpf(x) for x in a70["fit"]["coefficients"]]
c110 = [mp.mpf(x) for x in a110["fit"]["coefficients"]]
rep = {}
worst = mp.mpf(0)
for j, (x, y) in enumerate(zip(c70, c110)):
    diff = abs(x - y)
    rel = diff / max(1, abs(y))
    worst = max(worst, rel)
    print(f"c_{j}: |c70 - c110| = {mp.nstr(diff, 5)}   relative {mp.nstr(rel, 5)}   (dps {mp.mp.dps})")
    rep[f"c_{j}"] = {"abs_diff": mp.nstr(diff, 5), "rel_diff": mp.nstr(rel, 5)}
print("worst relative difference 70-run vs 110-run coefficients:", mp.nstr(worst, 5))
print("-log10(worst) =", mp.nstr(-mp.log10(worst), 6))

# last-fit vs first-fit inside the 110 run
cl = [mp.mpf(x) for x in a110["fit"]["coefficients_last_fit"]]
w2 = max(abs(x - y) / max(1, abs(y)) for x, y in zip(cl, c110))
print("110-run: worst relative difference first-fit vs last-fit:", mp.nstr(w2, 5))

# ---------- exploratory PSLQ for c_4 (leading coefficient) ----------
mp.mp.dps = 100
lam = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi ** 2
Li4 = mp.polylog(4, mp.mpf(1) / 2); Li5 = mp.polylog(5, mp.mpf(1) / 2); Li6 = mp.polylog(6, mp.mpf(1) / 2)
c4 = mp.mpf(a110["fit"]["coefficients"][4])
dfct = mp.mpf(a110["fit"]["k5_actual"]) if False else None
print("\nPSLQ target c_4 =", mp.nstr(c4, 60), " (dps", mp.mp.dps, ")")

basisA = [("lam^6", lam ** 6), ("pi^2 lam^4", pi2 * lam ** 4), ("pi^4 lam^2", pi2 ** 2 * lam ** 2), ("pi^6", pi2 ** 3),
          ("z3 lam^3", z3 * lam ** 3), ("pi^2 z3 lam", pi2 * z3 * lam), ("z3^2", z3 ** 2), ("z5 lam", z5 * lam)]
basisB = basisA + [("Li4(1/2) lam^2", Li4 * lam ** 2), ("pi^2 Li4(1/2)", pi2 * Li4),
                   ("Li5(1/2) lam", Li5 * lam), ("Li6(1/2)", Li6)]

def run(name, basis, target, label):
    vec = [target] + [v for _, v in basis]
    rel = mp.pslq(vec, maxcoeff=10 ** 8, maxsteps=10 ** 6)
    if rel is None:
        print(f"[{name}] {label}: no relation found (maxcoeff 1e8)")
        return None
    if rel[0] == 0:
        print(f"[{name}] {label}: relation without the target (basis dependent): {rel}")
        return None
    res = sum(r * v for r, v in zip(rel, vec))
    print(f"[{name}] {label}: relation {rel}")
    print(f"      residual of the relation = {mp.nstr(res, 5)} (dps {mp.mp.dps}); coefficients of the target:")
    for (nm, _), r in zip(basis, rel[1:]):
        if r != 0:
            print(f"        {nm}: {-mp.mpf(r) / rel[0]}")
    return rel

relA = run("A", basisA, c4, "c_4 vs {lam,pi^2,z3,z5} monomials")
relB = run("B", basisB, c4, "c_4 vs A + Li_n(1/2) monomials") if relA is None else None
