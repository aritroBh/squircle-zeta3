import sys
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from flint import fmpz_mat, fmpz
from math import gcd
from functools import reduce
from hpl_series import H_series, Li31, Li22

def find_relation_rows(vals, dps, margin=25, delta=0.99):
    mp.mp.dps = dps
    C = mp.mpf(10) ** (dps - margin)
    n = len(vals)
    M = fmpz_mat(n, n + 1)
    for i, v in enumerate(vals):
        M[i, i] = fmpz(1)
        M[i, n] = fmpz(int(mp.nint(C * v)))
    R = M.lll(delta=delta, eta=0.51)
    return [[R[r, c] for c in range(n + 1)] for r in range(n)]

def normalize(row):
    coeffs = [int(c) for c in row[:-1]]
    g = reduce(gcd, [abs(c) for c in coeffs if c])
    coeffs = [c // g for c in coeffs]
    if next((c for c in coeffs if c), 1) < 0: coeffs = [-c for c in coeffs]
    return coeffs

def try_rel(target, basis, names, dps, tag):
    rows = find_relation_rows([target] + basis, dps)
    coeffs = normalize(rows[0])
    resid = sum(c*v for c, v in zip(coeffs, [target] + basis))
    h = max(abs(c) for c in coeffs)
    ok = abs(resid) < mp.mpf(10)**(-(dps-60))
    print(f"{tag}: dps {dps} height {h} residual {mp.nstr(abs(resid),3)} {'RELATION' if ok else 'no'}")
    if ok:
        from fractions import Fraction as Fr
        a0 = coeffs[0]
        for i in range(1, len(coeffs)):
            if coeffs[i]:
                print(f"     {Fr(-coeffs[i], a0)!s:>8} * {names[i-1]}")
    return ok

mp.mp.dps = 500
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
w2_basis = [mp.pi**2, Li(2, mp.mpf(1)/3), l2**2, l2*l3, l3**2]
w2_names = ["pi2", "Li2(1/3)", "log2^2", "log2*log3", "log3^2"]
try_rel(Li(2, -mp.mpf(1)/2), w2_basis, w2_names, 400, "Li2(-1/2)")
try_rel(Li(2, mp.mpf(1)/4), w2_basis, w2_names, 400, "Li2(1/4)")
w3_basis = [mp.zeta(3), Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2), mp.pi**2*l2, mp.pi**2*l3,
            Li(2, mp.mpf(1)/3)*l2, Li(2, mp.mpf(1)/3)*l3, l2**3, l2**2*l3, l2*l3**2, l3**3]
w3_names = ["zeta3", "Li3(1/3)", "Li3(-1/2)", "pi2*log2", "pi2*log3",
            "Li2(1/3)*log2", "Li2(1/3)*log3", "log2^3", "log2^2*log3", "log2*log3^2", "log3^3"]
try_rel(Li(3, mp.mpf(2)/3), w3_basis, w3_names, 400, "Li3(2/3)")
try_rel(Li(3, mp.mpf(3)/4), w3_basis, w3_names, 400, "Li3(3/4)")
try_rel(Li(3, mp.mpf(1)/4), w3_basis, w3_names, 400, "Li3(1/4)")
