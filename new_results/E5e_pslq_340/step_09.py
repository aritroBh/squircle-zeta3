import json
from math import gcd
from functools import reduce
from fractions import Fraction
import mpmath as mp

C = json.load(open("consts360.json"))
rel = [7680000, 0, 1175040000, 49152000, 2574432000, -501811200, -111744000, 927640800,
       6229440, 397036800, 45175200, 287114250, 381542400, 77518080, 52715520, -1240896,
       -1751040, -10005440, -2977497, -24330240, -56401920]
labels = ["X"] + C["labels"]
strs = [C["vals"]["X"]] + [C["vals"][l] for l in C["labels"]]
maxcoeff = max(abs(c) for c in rel)
print("maxcoeff =", maxcoeff)

# ---- 340-digit residual ----
mp.mp.dps = 340
v = [mp.mpf(s) for s in strs]
res = sum(mp.mpf(c)*vv for c, vv in zip(rel, v))
print("residual @340 dps:", mp.nstr(res, 6))
bound = mp.mpf(maxcoeff) * mp.mpf(10)**(-330)
print("acceptance bound maxcoeff*1e-330 =", mp.nstr(bound, 6))
print("PASS:", abs(res) <= bound)

# residual with e2 instead of e1 (independent evaluation)
d = json.load(open('../E5d_nonlocal_pslq/b0_p5_convB_350.json'))
mp.mp.dps = 360
X2 = mp.mpf(d['e2']['P5'])/mp.sqrt(5)
mp.mp.dps = 340
v2 = [mp.mpf(mp.nstr(X2,350))] + v[1:]
res2 = sum(mp.mpf(c)*vv for c, vv in zip(rel, v2))
print("residual @340 with e2 leg:", mp.nstr(res2, 6), " PASS:", abs(res2) <= bound)

# ---- reduce ----
g = reduce(gcd, [abs(c) for c in rel if c != 0])
print("\ngcd of all coeffs:", g)
red = [c//g for c in rel]
print("reduced relation:", red)
print("reduced maxcoeff:", max(abs(c) for c in red))

# exact rational coordinates of X in this basis (X = sum_i q_i b_i, q_i = -rel[i+1]/rel[0])
q = [Fraction(-rel[i+1], rel[0]) for i in range(20)]
print("\nCoordinates of X = P/√5 in the 20-element basis:")
for lab, qq in zip(C["labels"], q):
    print(f"  {lab:16s} {qq}  = {float(qq):.10f}")
json.dump({"rel": rel, "red": red, "gcd": g,
           "q": [str(x) for x in q]}, open("relation.json","w"))
