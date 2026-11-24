"""w6 identification: coefficients c_0..c_4 from the dps-900 run (k=6..11): exact rational interpolation through k=6..10, check at k=11.
Writes coef_900.json (same format as coef_DPS.json).  All printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp
DPS = 900
mp.mp.dps = DPS
print("mp.mp.dps =", mp.mp.dps)
strs = {}
for k in range(6, 12):
    r = json.load(open(os.path.join(W, f"vals_{DPS}", f"k{k}_d5.json")))
    assert r["dps"] == DPS and r["k"] == k and r["d"] == 5
    strs[k] = r["value"]
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
def tomp(f): return mp.mpf(f.numerator) / mp.mpf(f.denominator)
Y = {k: Fr(strs[k]) for k in strs}
c1 = poly_through([(Fr(k), Y[k]) for k in range(6, 11)])
print("check residual of the quartic through k=6..10 at k=11:", mp.nstr(tomp(Y[11] - sum(c1[j] * Fr(11) ** j for j in range(5))), 5))
json.dump({"dps": DPS, "coefficients": [mp.nstr(tomp(c), DPS - 20) for c in c1]}, open(os.path.join(W, f"coef_{DPS}.json"), "w"), indent=1)
for j in range(5): print(f"c_{j} = {mp.nstr(tomp(c1[j]), 60)} ...")
