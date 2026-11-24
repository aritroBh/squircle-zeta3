# Step 4: EXACT (rational-arithmetic) structure of the final pieces as functions of k.
# For each key (i,j), pole r and Q(sqrt2)-component, the coefficient of 1/(v-r) is a sequence in k.
# We test exactly (Fractions) whether it is a polynomial in k, and for which k-range.
import sys, os, pickle, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import fastred
from fastred import QSI, ZERO, AQ, BQ

KS = list(range(2, 12))
fin = {k: pickle.load(open(f"final_{k}.pkl", "rb")) for k in KS}
name = {(ZERO.p0, ZERO.p1): "0", (AQ.p0, AQ.p1): "a", ((-AQ).p0, (-AQ).p1): "-a", (BQ.p0, BQ.p1): "b", ((-BQ).p0, (-BQ).p1): "-b"}

def coef(k, key, pole, comp):
    P = fin[k].get(key)
    if P is None: return Fr(0)
    for r, c in P[0].items():
        if name[(r.p0, r.p1)] == pole:
            return getattr(c, comp)
    return Fr(0)

def fdiff(seq, m):
    s = list(seq)
    for _ in range(m): s = [s[i+1]-s[i] for i in range(len(s)-1)]
    return s

def solve_poly(ks, vals):
    """exact Lagrange: coefficients c_0..c_{n-1} of the polynomial through (k, val)."""
    n = len(ks)
    A = [[Fr(k)**m for m in range(n)] + [Fr(v)] for k, v in zip(ks, vals)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        A[c] = [x / A[c][c] for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f*y for x, y in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]

def polyval(cs, k): return sum(c * Fr(k)**m for m, c in enumerate(cs))

summary = {}
for d in (1, 2, 3, 4):
    keys = sorted(set(kk for k in KS for kk in fin[k] if sum(kk) == d))
    degp = d - 1   # expected polynomial degree in k
    # find the smallest k0 such that, for all (key,pole,comp), the sequence on k0..11 is a polynomial of degree <= degp
    # (test: (degp+1)-th exact finite difference vanishes on every window inside k0..11)
    def poly_ok(k0, deg):
        ks = [k for k in KS if k >= k0]
        for key in keys:
            for pole in ("0", "a", "-a", "b", "-b"):
                for comp in ("p0", "p1", "p2", "p3"):
                    seq = [coef(k, key, pole, comp) for k in ks]
                    if any(fdiff(seq, deg+1)): return False
        return True
    k0 = next((k for k in KS if poly_ok(k, degp)), None)
    # also minimal degree on k >= k0+? (report minimal degree for the window k>=5 for d=4 etc.)
    mindeg = {}
    for k0_ in (2, 3, 4, 5):
        mindeg[k0_] = next((dg for dg in range(0, 8) if poly_ok(k0_, dg)), None)
    summary[d] = dict(weight=d+1, expected_degree=degp, first_k_with_exact_poly_of_that_degree=k0, min_exact_degree_for_k_ge=mindeg)
    print(f"weight {d+1} (key degree {d}): exact kernel coefficients polynomial of degree {degp} in k for all k >= {k0};"
          f" minimal exact degree on windows k>=2,3,4,5: {mindeg}", flush=True)
json.dump(summary, open("exact_poly_summary.json", "w"), indent=1)

# ---- weight 5 detail: cubic kernel through k=5..8, verified on k=9,10,11, defect at k=4 ----
d = 4
keys = sorted(set(kk for k in KS for kk in fin[k] if sum(kk) == d))
fitk = [5, 6, 7, 8]; chk = [9, 10, 11]
cub = {}   # (key,pole,comp) -> [c0..c3]
for key in keys:
    for pole in ("0", "a", "-a", "b", "-b"):
        for comp in ("p0", "p1"):
            cs = solve_poly(fitk, [coef(k, key, pole, comp) for k in fitk])
            for k in chk + [2, 3, 4]:
                pass
            ok = all(polyval(cs, k) == coef(k, key, pole, comp) for k in chk)
            assert ok, (key, pole, comp)
            cub[(key, pole, comp)] = cs
print("weight-5 kernels: cubic through k=5..8 reproduces k=9,10,11 EXACTLY for every (key,pole,component)")
defect = {}
for kk_, cs in cub.items():
    dv = coef(4, kk_[0], kk_[1], kk_[2]) - polyval(cs, 4)
    if dv != 0: defect[kk_] = dv
print("k=4 defect kernel (actual K^(4) minus cubic-extrapolated), coefficients of 1/(v-pole) in comp p0 + p1*sqrt2:")
for kk_, dv in sorted(defect.items()): print("   ", kk_, dv)
print("cubic kernel coefficients c_m (m=0..3) per (key,pole,comp):")
for kk_, cs in sorted(cub.items()):
    if any(cs): print("   ", kk_, [str(c) for c in cs])
pickle.dump(dict(cub=cub, defect=defect, keys=keys), open("exact_cubic_kernel_w5.pkl", "wb"))
