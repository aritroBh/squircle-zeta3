"""Fixed-basis LLL for W4'' : the SAME 60 numbers as the W4 run (m3_pslq_1150.out kept Li4's).
Usage: step_pslq2.py <W4pp file> <D>."""
import os, sys, time, itertools
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
import flint
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B
import polylog_exact as PE

WFILE, D = sys.argv[1], int(sys.argv[2])
mp.mp.dps = D + 40
T0 = time.time()
W = mp.mpf(open(WFILE).read().strip())

r2 = mp.sqrt(2); lam = mp.log(2); L = mp.log(1+r2); b = r2-1
pi2 = mp.pi**2; Li2b = mp.polylog(2, b)
z3 = mp.zeta(3); Li3b = mp.polylog(3, b); Li3mb = mp.polylog(3, -b); Li3h = mp.polylog(3, 1/r2)
w1 = [("log2", lam), ("L", L)]
w2 = [("pi^2", pi2), ("Li2(b)", Li2b)]
w3 = [("z3", z3), ("Li3(b)", Li3b), ("Li3(-b)", Li3mb), ("Li3(1/sqrt2)", Li3h)]
prods = []
for c in itertools.combinations_with_replacement(range(2), 4):
    prods.append(("*".join(w1[i][0] for i in c), mp.fprod(w1[i][1] for i in c)))
for (n2, v2) in w2:
    for c in itertools.combinations_with_replacement(range(2), 2):
        prods.append((n2+"*"+"*".join(w1[i][0] for i in c), v2*mp.fprod(w1[i][1] for i in c)))
for c in itertools.combinations_with_replacement(range(2), 2):
    prods.append(("*".join(w2[i][0] for i in c), mp.fprod(w2[i][1] for i in c)))
for (n3, v3) in w3:
    for (n1, v1) in w1:
        prods.append((n3+"*"+n1, v3*v1))

KEPT_Z = [QS2(Fr(1, 2)), B, -B, QS2(0, Fr(1, 2)), -B*B, QS2(2)-SQ2,
          QS2(1)-QS2(0, Fr(1, 2)), QS2(-1, Fr(1, 2))]
kept = [(f"Li4({z})", mp.polylog(4, PE.DE.float_exact(z))) for z in KEPT_Z]

def qs2(lst):
    out = []
    for nm, v in lst:
        out += [(nm, v), ("sqrt2*"+nm, r2*v)]
    return out

def lll(vals, digits):
    C = mp.mpf(10)**digits
    n = len(vals)
    rows = [[1 if i == j else 0 for j in range(n)] + [int(mp.nint(C*vals[i]))] for i in range(n)]
    Bm = flint.fmpz_mat(rows).lll()
    best = [int(Bm[0, j]) for j in range(n)]
    with mp.workdps(digits + 20):
        resid = mp.fsum(c*v for c, v in zip(best, vals))
    return best, max(abs(c) for c in best), resid

full = [("W4''", W)] + qs2(prods) + qs2(kept)
print(f"basis size = {len(full)} (should be 61)")
vec, mx, res = lll([v for _, v in full], D)
thresh = 10**(D/(2*len(full)))
print(f"LLL at {D} digits: W4'' coef = {vec[0]}, max|coef| = {mx:.3e}, residual = {mp.nstr(res,3)}")
print(f"  (random size ~ 1e{D/len(full):.1f}; threshold {thresh:.1e})")
import json
json.dump({"D": D, "vec0": vec[0], "maxcoef": str(mx), "residual": mp.nstr(res, 6),
           "threshold": f"{thresh:.2e}", "kept": [n for n, _ in kept],
           "relation": {nm: c for c, nm in zip(vec, [n for n, _ in full]) if c}},
          open(f"pslq2_w4pp_{D}.json", "w"), indent=1)
rel = " ".join(f"{c:+d}*{nm}" for c, nm in zip(vec, [n for n, _ in full]) if c)
print("RELATION:", rel)
print(f"saved pslq2_w4pp_{D}.json; {time.time()-T0:.0f}s")
