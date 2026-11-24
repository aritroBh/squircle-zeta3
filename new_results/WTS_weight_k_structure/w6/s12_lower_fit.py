"""w6 identification step 5 (structure guide): polynomial laws for W_4(k) (d=3) and W_5(k) (d=4) from vals_150
(dps 150, k=4..15).  Exact rational interpolation of the decimal strings; smallest degree from finite differences
on k >= weight; coefficients saved to lower_coef.json.  All printed numbers come from this script."""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp
DPS = 150
mp.mp.dps = DPS
print("mp.mp.dps =", mp.mp.dps)
def tomp(f): return mp.mpf(f.numerator) / mp.mpf(f.denominator)
def diffs(seq, n):
    cur = list(seq)
    for _ in range(n): cur = [cur[i + 1] - cur[i] for i in range(len(cur) - 1)]
    return cur
def poly_through(pts):
    n = len(pts); coeffs = [Fr(0)] * n
    for a in range(n):
        ka, ya = pts[a]; num = [Fr(1)]; den = Fr(1)
        for c in range(n):
            if c == a: continue
            kc = pts[c][0]; new = [Fr(0)] * (len(num) + 1)
            for t, v in enumerate(num): new[t + 1] += v; new[t] -= kc * v
            num = new; den *= (ka - kc)
        for t, v in enumerate(num): coeffs[t] += ya * v / den
    return coeffs
out = {}
for d in (3, 4):
    wt = d + 1
    Y = {k: Fr(json.load(open(os.path.join(W, "vals_150", f"k{k}_d{d}.json")))["value"]) for k in range(4, 16)}
    for kmin in sorted(set((max(4, wt - 1), wt))):
        ks = list(range(kmin, 16))
        ys = [Y[k] for k in ks]
        deg = None
        for n in range(1, len(ys)):
            m = max(abs(tomp(x)) for x in diffs(ys, n))
            if m < mp.mpf(10) ** -(DPS - 25):
                deg = n - 1; mres = m; break
        print(f"weight {wt}, window k={kmin}..15: smallest degree with vanishing (d+1)-th difference: {deg}" + (f" (residual {mp.nstr(mres,4)})" if deg is not None else " (none)"))
        if deg is not None and kmin == wt:
            c = poly_through([(Fr(k), Y[k]) for k in ks[:deg + 1]])
            chk = max(abs(tomp(Y[k] - sum(c[j] * Fr(k) ** j for j in range(deg + 1)))) for k in ks)
            print(f"   fit through k={ks[:deg+1]}, max residual over window = {mp.nstr(chk,4)}")
            out[f"W{wt}"] = [mp.nstr(tomp(x), DPS - 10) for x in c]
            for j, x in enumerate(c): print(f"   c_{j} = {mp.nstr(tomp(x), 60)}")
            kb = wt - 1
            pk = sum(c[j] * Fr(kb) ** j for j in range(deg + 1))
            if kb in Y: print(f"   k={kb}: W - poly = {mp.nstr(tomp(Y[kb] - pk), 8)}")
json.dump(out, open(os.path.join(W, "lower_coef.json"), "w"), indent=1)
