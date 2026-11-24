import subprocess, sys
code = r'''"""Task 5: m3_pslq-style LLL for W4'' with the same 60-number basis as W4.
Usage: step_pslq.py <W4pp value file> <D>."""
import os, sys, time, itertools
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
import flint
from fractions import Fraction as Fr
from qsqrt2 import QS2, SQ2, A, B, sunit_log
import polylog_exact as PE

WFILE, D = sys.argv[1], int(sys.argv[2])
mp.mp.dps = D + 40
T0 = time.time()
W = mp.mpf(open(WFILE).read().strip())
print(f"W4'' from {WFILE}; LLL at {D} digits")

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
print(f"{len(prods)} weight-4 products ({time.time()-T0:.0f}s)", flush=True)

cands = set()
for j in range(-6, 7):
    for m in range(-4, 5):
        for sg in (1, -1):
            z = QS2(sg)*SQ2**j*A**m
            if not (z < QS2(1)) or z == QS2(0):
                continue
            try:
                sunit_log(QS2(1)-z)
            except ValueError:
                continue
            if abs(float(z)) > 1:
                continue
            cands.add(z)
cands = sorted(cands, key=lambda z: (abs(float(z)) > 0.9, -abs(float(z))))
print(f"{len(cands)} candidate Li4 arguments ({time.time()-T0:.0f}s)", flush=True)

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

base = qs2(prods)
kept = []
for z in cands:
    val = mp.polylog(4, PE.DE.float_exact(z))
    trial = base + qs2(kept) + qs2([(f"Li4({z})", val)])
    DS = min(D, 1000)
    vec, mx, res = lll([v for _, v in trial], DS)
    thresh = 10**(DS/(2*len(trial)))
    dependent = mx < thresh and (vec[-1] != 0 or vec[-2] != 0)
    print(f"  Li4({z}): {'dependent' if dependent else 'kept'} (max coef {mx:.2e}) ({time.time()-T0:.0f}s)", flush=True)
    if not dependent:
        kept.append((f"Li4({z})", val))
print(f"kept {len(kept)} independent Li4's:", [n for n, _ in kept], flush=True)

full = [("W4''", W)] + qs2(prods) + qs2(kept)
vec, mx, res = lll([v for _, v in full], D)
thresh = 10**(D/(2*len(full)))
print(f"LLL on W4'' + {len(full)-1} basis numbers at {D} digits: max|coef| = {mx:.3e}, residual = {mp.nstr(res,3)}")
print(f"  (random size ~ 1e{D/len(full):.1f}; genuine-relation threshold {thresh:.1e})")
found = vec[0] != 0 and mx < thresh
print("W4'' coefficient:", vec[0])
if found:
    rel = " ".join(f"{c:+d}*{nm}" for c, nm in zip(vec, [n for n, _ in full]) if c)
    print("RELATION FOUND:")
    print(rel)
import json
out = {"D": D, "vec0": vec[0], "maxcoef": str(mx), "residual": mp.nstr(res, 6),
       "threshold": f"{thresh:.2e}", "found": bool(found),
       "kept": [n for n, _ in kept],
       "relation": {nm: c for c, nm in zip(vec, [n for n, _ in full]) if c}}
json.dump(out, open(f"pslq_w4pp_{D}.json", "w"), indent=1)
print(f"saved pslq_w4pp_{D}.json; total {time.time()-T0:.0f}s")
'''
open("step_pslq.py", "w").write(code)
r = subprocess.run([sys.executable, "step_pslq.py", "W4pp_800.txt", "700"], capture_output=True, text=True, timeout=2700)
print(r.stdout[-4500:])
print("ERR", r.stderr[-600:])
