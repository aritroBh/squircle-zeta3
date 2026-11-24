"""w6 identification step 15: exact-interpolation coefficients of W_4(k) (quadratic, k=4..8) and W_5(k) (cubic, k=5..10) from the dps-400 data.
Writes lower_coef_400.json.  All printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp
DPS = 400
mp.mp.dps = DPS
print("mp.mp.dps =", mp.mp.dps)
def tomp(f): return mp.mpf(f.numerator) / mp.mpf(f.denominator)
def poly_through(pts):
    n = len(pts); co = [Fr(0)] * n
    for a in range(n):
        ka, ya = pts[a]; num = [Fr(1)]; den = Fr(1)
        for c in range(n):
            if c == a: continue
            kc = pts[c][0]; new = [Fr(0)] * (len(num) + 1)
            for t, v in enumerate(num): new[t + 1] += v; new[t] -= kc * v
            num = new; den *= (ka - kc)
        for t, v in enumerate(num): co[t] += ya * v / den
    return co
out = {}
for d, ks, deg in ((3, list(range(4, 9)), 2), (4, list(range(5, 11)), 3)):
    Y = {k: Fr(json.load(open(os.path.join(W, "vals_400", f"k{k}_d{d}.json")))["value"]) for k in ks}
    c = poly_through([(Fr(k), Y[k]) for k in ks[:deg + 1]])
    chk = [tomp(Y[k] - sum(c[j] * Fr(k) ** j for j in range(deg + 1))) for k in ks[deg + 1:]]
    print(f"weight {d+1}: degree {deg} fit through k={ks[:deg+1]}; check residuals at k={ks[deg+1:]}: {[mp.nstr(x, 3) for x in chk]}")
    out[f"W{d+1}"] = [mp.nstr(tomp(x), DPS - 15) for x in c]
json.dump(out, open(os.path.join(W, "lower_coef_400.json"), "w"), indent=1)
