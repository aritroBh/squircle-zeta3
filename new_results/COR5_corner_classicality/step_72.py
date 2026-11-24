import sys
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from flint import fmpz_mat, fmpz
from math import gcd
from functools import reduce

def find_relation_rows(vals, dps, margin=30, delta=0.99):
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
    g = reduce(gcd, [abs(c) for c in coeffs if c]) or 1
    coeffs = [c // g for c in coeffs]
    if next((c for c in coeffs if c), 1) < 0: coeffs = [-c for c in coeffs]
    return coeffs

mp.mp.dps = 1120
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
l313, l31m = Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2)
z3 = mp.zeta(3); l213 = Li(2, mp.mpf(1)/3); pi2 = mp.pi**2
def mk_basis():
    return [Li(4, mp.mpf(1)/3), Li(4, mp.mpf(2)/3), Li(4, mp.mpf(3)/4), Li(4, mp.mpf(1)/2), Li(4, mp.mpf(1)/9),
            l313*l2, l313*l3, l31m*l2, l31m*l3, z3*l2, z3*l3,
            l213**2, l213*l2**2, l213*l2*l3, l213*l3**2, pi2*l213,
            pi2*l2**2, pi2*l2*l3, pi2*l3**2, mp.pi**4,
            l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4]
B = mk_basis()
t = Li(4, -mp.mpf(1)/2)
for dps in (700, 1000):
    rows = find_relation_rows([t] + B, dps)
    coeffs = normalize(rows[0])
    mp.mp.dps = 1120
    resid = sum(c*v for c, v in zip(coeffs, [t] + B))
    print(f"Li4(-1/2) dps {dps}: height {max(abs(c) for c in coeffs)}, resid@1120 {mp.nstr(abs(resid),3)}")
