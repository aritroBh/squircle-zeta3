import mpmath as mp
from flint import fmpz_mat, fmpz
from fractions import Fraction as Fr
from math import gcd
from functools import reduce

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
    if next(c for c in coeffs if c) < 0:
        coeffs = [-c for c in coeffs]
    return coeffs

mp.mp.dps = 1050
W3 = mp.mpf(open("hex_W3.txt").read())
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
basis3 = [mp.zeta(3), Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2), mp.pi**2*l2, mp.pi**2*l3,
          Li(2, mp.mpf(1)/3)*l2, Li(2, mp.mpf(1)/3)*l3, l2**3, l2**2*l3, l2*l3**2, l3**3]
cand = {}
for dps in (600, 900):
    rows = find_relation_rows([W3] + basis3, dps)
    coeffs = normalize(rows[0])
    resid = sum(c*v for c, v in zip(coeffs, [W3] + basis3))
    print(f"dps {dps}: normalized height {max(abs(c) for c in coeffs)}, residual {mp.nstr(abs(resid), 3)}")
    cand[dps] = coeffs
print("IDENTICAL normalized relation:", cand[600] == cand[900])
