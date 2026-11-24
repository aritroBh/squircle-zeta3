"""w6 identification step 1: coefficients c_0..c_4 of W_6(k) (k>=6) from the high-precision reruns.
usage: s08_fit_hi.py DPS      (reads vals_DPS/k{6..15}_d5.json written by s02_eval.py --dps DPS)
Exact rational interpolation of the decimal value strings through k=6..10 ("first fit") and k=11..15 ("last fit"),
check residuals of the first fit at k=11..15, difference first-vs-last fit.  Writes coef_DPS.json.
All printed numbers come from this script.
"""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp

DPS = int(sys.argv[1])
mp.mp.dps = DPS          # AFTER all imports
print("mp.mp.dps =", mp.mp.dps)
strs = {}
for k in range(6, 16):
    r = json.load(open(os.path.join(W, f"vals_{DPS}", f"k{k}_d5.json")))
    assert r["dps"] == DPS and r["k"] == k and r["d"] == 5
    strs[k] = r["value"]

def poly_through(pts):
    n = len(pts); coeffs = [Fr(0)] * n
    for a in range(n):
        ka, ya = pts[a]; num = [Fr(1)]; den = Fr(1)
        for c in range(n):
            if c == a: continue
            kc = pts[c][0]
            new = [Fr(0)] * (len(num) + 1)
            for t, v in enumerate(num):
                new[t + 1] += v; new[t] -= kc * v
            num = new; den *= (ka - kc)
        for t, v in enumerate(num):
            coeffs[t] += ya * v / den
    return coeffs

def tomp(f): return mp.mpf(f.numerator) / mp.mpf(f.denominator)
Y = {k: Fr(strs[k]) for k in strs}
c1 = poly_through([(Fr(k), Y[k]) for k in range(6, 11)])
c2 = poly_through([(Fr(k), Y[k]) for k in range(11, 16)])
p1 = lambda k: sum(c1[j] * Fr(k) ** j for j in range(5))
print("check residuals of first fit (k=6..10) at k=11..15:")
for k in range(11, 16):
    print(f"   k={k}: W_6(k) - p(k) = {mp.nstr(tomp(Y[k] - p1(k)), 5)}")
rel = max(abs(tomp(a) - tomp(b)) / max(1, abs(tomp(a))) for a, b in zip(c1, c2))
print("worst relative |first-fit - last-fit| =", mp.nstr(rel, 5))
for j in range(5):
    print(f"c_{j} = {mp.nstr(tomp(c1[j]), DPS - 12)}")
json.dump({"dps": DPS, "coefficients": [mp.nstr(tomp(c), DPS - 12) for c in c1],
           "coefficients_last_fit": [mp.nstr(tomp(c), DPS - 12) for c in c2],
           "worst_rel_first_vs_last": mp.nstr(rel, 5)}, open(os.path.join(W, f"coef_{DPS}.json"), "w"), indent=1)
