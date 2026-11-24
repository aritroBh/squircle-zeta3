"""Assemble values_w4.json from the saved quadrature results (dps 90 and dps 120) and the exact-fit checks."""
import json
from w4_common import *
from w4_known import *
set_dps(120); print("mp.mp.dps =", mp.mp.dps)
def val(w, k, dp): return mp.mpf(json.load(open(os.path.join(RES, f"w{w}_k{k}_dps{dp}.json")))["value"])
ks = list(range(3, 11))
v120 = {k: val(4, k, 120) for k in ks}; v90 = {k: val(4, k, 90) for k in ks}
print("max |W_4(k)[dps90] - W_4(k)[dps120]| over k=3..10:", mp.nstr(max(abs(v90[k] - v120[k]) for k in ks), 5))
def fit(kk):
    A = mp.matrix(3, 3); y = mp.matrix(3, 1)
    for r, k in enumerate(kk):
        for c in range(3): A[r, c] = mp.mpf(k)**c
        y[r] = v120[k]
    return mp.lu_solve(A, y)
cA = fit([4, 5, 6]); cB = fit([8, 9, 10]); cC = fit([5, 7, 9])
print("fit(4,5,6) vs fit(8,9,10) max coefficient difference:", mp.nstr(max(abs(cA[j] - cB[j]) for j in range(3)), 5))
print("fit(4,5,6) vs fit(5,7,9)  max coefficient difference:", mp.nstr(max(abs(cA[j] - cC[j]) for j in range(3)), 5))
P = lambda k: sum(cA[j] * mp.mpf(k)**j for j in range(3))
resid = {k: P(k) - v120[k] for k in ks}
for k in ks: print(f"  k={k}: P(k) - W_4(k) = {mp.nstr(resid[k], 6)}")
# third differences
d3 = [v120[k+3] - 3*v120[k+2] + 3*v120[k+1] - v120[k] for k in range(3, 8)]
print("third differences starting at k=3..7:", [mp.nstr(x, 6) for x in d3])
# closed form of c_2 from PSLQ (exploratory) re-verified at both precisions
lam = mp.log(2); z3 = mp.zeta(3)
c2_cf = -(3*mp.pi**4 + 20*mp.pi**2*lam**2 + 20*lam**4 + 120*z3*lam) / 960
print("c_2 - closed form (dps 120 fit):", mp.nstr(cA[2] - c2_cf, 5))
# same check from the dps-90 data
set_dps(90)
A = mp.matrix(3, 3); y = mp.matrix(3, 1)
for r, k in enumerate([4, 5, 6]):
    for c in range(3): A[r, c] = mp.mpf(k)**c
    y[r] = v90[k]
c90 = mp.lu_solve(A, y)
lam = mp.log(2); z3 = mp.zeta(3)
print("c_2(dps90 fit) - closed form, mp.mp.dps =", mp.mp.dps, ":", mp.nstr(c90[2] - (-(3*mp.pi**4 + 20*mp.pi**2*lam**2 + 20*lam**4 + 120*z3*lam)/960), 5))
set_dps(120)
delta = P(3) - v120[3]
out = {
  "weight": 4, "polynomial_law": {"valid_for": "k >= 4", "degree": 2, "fit_points": [4, 5, 6],
         "coeffs_c0_c1_c2": [mp.nstr(cA[j], 112) for j in range(3)],
         "residual_at_k": {str(k): mp.nstr(resid[k], 6) for k in ks},
         "fit_456_vs_8910_max_coeff_diff": mp.nstr(max(abs(cA[j] - cB[j]) for j in range(3)), 5)},
  "k3_exception": {"P(3)-W_4(3)": mp.nstr(delta, 112), "third_diff_from_k3": mp.nstr(d3[0], 112)},
  "c2_closed_form": "-(3 pi^4 + 20 pi^2 log^2 2 + 20 log^4 2 + 120 zeta(3) log 2)/960   [PSLQ, exploratory]",
  "values_dps120": {str(k): mp.nstr(v120[k], 112) for k in ks},
  "values_dps90": {str(k): mp.nstr(v90[k], 85) for k in ks},
  "W_4(2)": "0 (no weight-4 piece in M_2: top weight of M_2 is 3)",
  "controls_dps90": json.load(open(os.path.join(RES, "controls_dps90.json"))),
  "controls_dps120": json.load(open(os.path.join(RES, "controls_dps120.json"))),
}
json.dump(out, open(os.path.join(HERE, "values_w4.json"), "w"), indent=1)
print("wrote values_w4.json")
