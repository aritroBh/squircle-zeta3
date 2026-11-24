"""Exact finite-difference/polynomial fit of W_4(k) = c0 + c1 k + c2 k^2 (k>=4) at a given dps, plus controls.
Usage: nice -n 10 python id_fit.py <dps>.  Reads results/w4_k<k>_dps<dps>.json (k=3..8).  Writes results/id_coeffs_dps<dps>.json.
PRECISION: mp.mp.dps is set AFTER all imports (w4_known imports nothing that changes it, but we set it anyway)."""
import sys, json
from w4_common import *
from w4_known import *
dps = int(sys.argv[1]); set_dps(dps); print("mp.mp.dps =", mp.mp.dps)
def val(k): return mp.mpf(json.load(open(os.path.join(RES, f"w4_k{k}_dps{dps}.json")))["value"])
v = {k: val(k) for k in range(3, 9)}
# controls
w43 = W4_thm(); w44 = W4pp_txt() if dps <= 790 else None
print("control W_4(3) - thm:W4          =", mp.nstr(v[3] - w43, 5), " (mp.mp.dps", mp.mp.dps, ")")
print("control W_4(4) - W4pp_800.txt    =", mp.nstr(v[4] - w44, 5), " (txt has 800 digits)")
print("control W_4(4) - w4pp closed form=", mp.nstr(v[4] - W4pp_closed(), 5))
def fit(kk):
    A = mp.matrix(3, 3); y = mp.matrix(3, 1)
    for r, k in enumerate(kk):
        for c in range(3): A[r, c] = mp.mpf(k)**c
        y[r] = v[k]
    return mp.lu_solve(A, y)
cA = fit([4, 5, 6]); cB = fit([6, 7, 8]); cC = fit([4, 6, 8])
print("max coeff diff fit(4,5,6) vs fit(6,7,8):", mp.nstr(max(abs(cA[j] - cB[j]) for j in range(3)), 5))
print("max coeff diff fit(4,5,6) vs fit(4,6,8):", mp.nstr(max(abs(cA[j] - cC[j]) for j in range(3)), 5))
P = lambda k: sum(cA[j] * mp.mpf(k)**j for j in range(3))
for k in (7, 8): print(f"  P({k}) - W_4({k}) =", mp.nstr(P(k) - v[k], 5))
delta = P(3) - v[3]
print("delta = P(3)-W_4(3) =", mp.nstr(delta, 40))
for j in range(3): print(f"c{j} =", mp.nstr(cA[j], 60))
json.dump({"dps": dps, "mp_dps": mp.mp.dps, "c": [mp.nstr(cA[j], dps - 5) for j in range(3)],
           "delta": mp.nstr(delta, dps - 5), "W4": {str(k): mp.nstr(v[k], dps - 5) for k in v},
           "fit_diff_456_vs_678": mp.nstr(max(abs(cA[j] - cB[j]) for j in range(3)), 5)},
          open(os.path.join(RES, f"id_coeffs_dps{dps}.json"), "w"), indent=1)
