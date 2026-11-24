"""Finite-difference / polynomial-law analysis.  Usage: python w4_poly.py <dps> <weight> <kmin> <kmax>
Reads results/w<weight>_k<k>_dps<dps>.json.  Prints the difference table, the smallest degree d with vanishing
(d+1)-th differences, and the exact-through-first-(d+1)-points fit with residuals at the remaining k."""
import sys, json
from w4_common import *
dps, w, kmin, kmax = (int(x) for x in sys.argv[1:5])
set_dps(dps)
print("mp.mp.dps =", mp.mp.dps)
ks = list(range(kmin, kmax + 1))
vals = [mp.mpf(json.load(open(os.path.join(RES, f"w{w}_k{k}_dps{dps}.json")))["value"]) for k in ks]
for k, v in zip(ks, vals): print(f"  W_{w}({k}) = {mp.nstr(v, 40)}")
# difference table
tab = [vals]
while len(tab[-1]) > 1:
    p = tab[-1]; tab.append([p[i + 1] - p[i] for i in range(len(p) - 1)])
print("max |n-th finite difference| (n = order):")
mx = []
for n, row in enumerate(tab):
    m = max(abs(x) for x in row); mx.append(m)
    print(f"  order {n}: n_terms {len(row)}  max|.| = {mp.nstr(m, 6)}")
tol = mp.mpf(10) ** (-(dps - 15))
d = -1
for n in range(1, len(tab)):
    if mx[n] < tol: d = n - 1; break
print("working-precision threshold:", mp.nstr(tol, 3))
print("smallest degree d with vanishing (d+1)-th differences:", d)
if d >= 0:
    # exact interpolation through first d+1 points (Newton form -> monomial coefficients via linear solve)
    A = mp.matrix(d + 1, d + 1); y = mp.matrix(d + 1, 1)
    for r in range(d + 1):
        for c in range(d + 1): A[r, c] = mp.mpf(ks[r]) ** c
        y[r] = vals[r]
    cj = mp.lu_solve(A, y)
    for j in range(d + 1): print(f"  c_{j} = {mp.nstr(cj[j], dps - 5)}")
    print("check at remaining k (fit - computed), mp.mp.dps =", mp.mp.dps)
    worst = mp.mpf(0)
    for k, v in zip(ks[d + 1:], vals[d + 1:]):
        fit = sum(cj[j] * mp.mpf(k) ** j for j in range(d + 1))
        r = fit - v; worst = max(worst, abs(r))
        print(f"    k={k}: residual {mp.nstr(r, 5)}")
    print("worst residual:", mp.nstr(worst, 5))
    json.dump({"weight": w, "dps": dps, "ks": ks, "degree": d, "coeffs": [mp.nstr(cj[j], dps - 5) for j in range(d + 1)],
               "worst_residual": mp.nstr(worst, 5), "diff_max": [mp.nstr(m, 6) for m in mx]},
              open(os.path.join(RES, f"poly_w{w}_dps{dps}_k{kmin}-{kmax}.json"), "w"), indent=1)
