"""verify step: coefficients c_0..c_4 from OWN 110-digit values (quartic through k=6..10, exact rational interpolation)
versus the stored fits coef_110/220/250/400/900.json (their runs).  usage: verify_coef_compare.py OWNDIR"""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import mpmath as mp
OWN = sys.argv[1]
mp.mp.dps = 130
print("mp.mp.dps =", mp.mp.dps)
Y = {k: Fr(json.load(open(os.path.join(OWN, f"own_k{k}_d5_110.json")))["value"]) for k in range(6, 13)}
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
tomp = lambda f: mp.mpf(f.numerator) / mp.mpf(f.denominator)
c = poly_through([(Fr(k), Y[k]) for k in range(6, 11)])
c2 = poly_through([(Fr(k), Y[k]) for k in range(8, 13)])
print("own fit k=6..10 vs own fit k=8..12, relative differences:", [mp.nstr(abs(tomp(a) - tomp(b)) / max(1, abs(tomp(a))), 3) for a, b in zip(c, c2)])
for name in ("analysis_vals_110.json", "coef_220.json", "coef_250.json", "coef_400.json", "coef_900.json"):
    d = json.load(open(os.path.join(W, name)))
    th = d["fit"]["coefficients"] if "fit" in d else d["coefficients"]
    print(name, "| (own - theirs) per coefficient:", [mp.nstr(tomp(c[j]) - mp.mpf(th[j]), 3) for j in range(5)])
for j in range(5):
    print(f"c_{j} (own) = {mp.nstr(tomp(c[j]), 105)}")
