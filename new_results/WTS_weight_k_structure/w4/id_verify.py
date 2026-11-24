"""Independent verification of the closed forms: evaluate c_2, c_1, c_0 from the exact coordinates in results/id_lll_normalised.json
(atoms recomputed here from scratch) and compare with (a) the coefficients obtained from the W_4(k) quadrature route (id_fit.py: results/id_coeffs_dps360.json),
(b) the predicted W_4(k) = c0 + c1 k + c2 k^2 against the quadrature values W_4(k), k = 4..8.   PRECISION: mp.mp.dps set after imports and printed."""
import json, sys
from fractions import Fraction as F
from id_lib import *
dps = int(sys.argv[1]) if len(sys.argv) > 1 else 360
set_dps(dps); print("mp.mp.dps =", mp.mp.dps)
norm = json.load(open(os.path.join(RES, "id_lll_normalised.json")))
S = with_sqrt2(atoms())
def cf(nm):
    rel = next(iter(norm[nm].values()))
    return sum(frac_mp(F(v)) * S[k] for k, v in rel.items())
ref = json.load(open(os.path.join(RES, f"id_coeffs_dps{dps}.json")))
cc = [cf("c0"), cf("c1"), cf("c2")]
for j in range(3):
    print(f"c{j}: closed form - quadrature-fit value = {mp.nstr(cc[j] - mp.mpf(ref['c'][j]), 5)}   (dps {mp.mp.dps})")
for k in range(4, 9):
    pred = cc[0] + cc[1] * k + cc[2] * k * k
    print(f"  k={k}: closed-form W_4(k) - quadrature W_4(k) = {mp.nstr(pred - mp.mpf(ref['W4'][str(k)]), 5)}")
print("c0 =", mp.nstr(cc[0], 40)); print("c1 =", mp.nstr(cc[1], 40)); print("c2 =", mp.nstr(cc[2], 40))
