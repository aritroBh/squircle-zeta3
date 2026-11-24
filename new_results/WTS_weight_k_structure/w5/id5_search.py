# Identification step 5: staged LLL searches (exact Gram) for c_0,c_1,c_2,delta and the 10 block integrals
# in families of weight-5 constants.  usage: id5_search.py <family> <D>
import sys, os, json, time, itertools
sys.dont_write_bytecode = True
from id_lib import *
fam = sys.argv[1]; D = int(sys.argv[2])
setprec(1450)
C = Consts()
R = json.load(open("results/cm_D1500_m10.json"))
print("flint ctx.prec =", ctx.prec, "bits; search digits D =", D, "family", fam, flush=True)
targets = {f"c{i}": arb(R["c"][i]) for i in range(4)}
targets["delta"] = arb(R["delta"])
B = json.load(open("results/blocks_D1500_m10.json"))["blocks"]
for k, v in B.items(): targets["blk_" + k] = arb(v["mid"])
lam, L, pi, z3, z5, Li4h, Li5h, sq2 = C.lam, C.L, C.pi, C.z3, C.z5, C.Li4h, C.Li5h, C.sq2

def monomials(gens, wt):
    """all monomials of total weight wt in gens=[(name,value,weight)]"""
    out = []
    def rec(i, rem, name, val):
        if rem == 0: out.append((name, val)); return
        if i == len(gens): return
        g, v, w = gens[i]
        e = 0; nm = name; vv = val
        while e * w <= rem:
            rec(i + 1, rem - e * w, nm, vv)
            e += 1; nm = nm + ("" if e == 0 else f"*{g}"); vv = vv * v
    rec(0, wt, "", arb(1))
    return [(n.lstrip("*") or "1", v) for n, v in out]

def family(name):
    if name == "F8":
        return [("z5", z5), ("pi2 z3", pi**2*z3), ("z3 lam2", z3*lam**2), ("pi4 lam", pi**4*lam), ("pi2 lam3", pi**2*lam**3), ("lam5", lam**5), ("Li4h lam", Li4h*lam), ("Li5h", Li5h)]
    if name == "F8s":
        f = family("F8"); return f + [("sq2*" + n, v*sq2) for n, v in f]
    gens20 = [("lam", lam, 1), ("L", L, 1), ("pi2", pi**2, 2), ("z3", z3, 3), ("Li4h", Li4h, 4), ("z5", z5, 5), ("Li5h", Li5h, 5)]
    if name == "F20":
        return monomials(gens20, 5)
    if name == "F20s":
        f = family("F20"); return f + [("sq2*" + n, v*sq2) for n, v in f]
    raise SystemExit("unknown family")

basis = family(fam)
print("basis size", len(basis), flush=True)
only = sys.argv[3:] 
res = {}
for tname, x in targets.items():
    if only and tname not in only: continue
    t0 = time.time()
    best, top = lll_relation([x] + [v for _, v in basis], D)
    coeffs, last, mx = best
    # the relation row of interest must have a nonzero coefficient on the target
    ok = coeffs[0] != 0 and abs(last) < 10**8 and mx < 10 ** (D // (len(basis) + 1) - 4)
    r = residual(coeffs, [x] + [v for _, v in basis])
    print(f"{tname:16s} n={len(basis)+1} maxcoef={mx:.3e} lastentry={last} resid~{r.mid().str(5)} -> {'RELATION' if ok else 'none'}  ({time.time()-t0:.1f}s)", flush=True)
    if ok:
        print("   coeffs:", {n: c for n, c in zip(["target"] + [n for n, _ in basis], coeffs) if c != 0}, flush=True)
    res[tname] = dict(coeffs=coeffs, last=last, mx=str(mx), ok=ok)
json.dump(res, open(f"results/id5_{fam}_D{D}.json", "w"), indent=1)
