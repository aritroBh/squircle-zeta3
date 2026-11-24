import mpmath as mp
from flint import fmpz_mat, fmpz
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import time

def find_relation_rows(vals, dps, margin=25, delta=0.99):
    mp.mp.dps = dps
    C = mp.mpf(10) ** (dps - margin)
    n = len(vals)
    M = fmpz_mat(n, n + 1)
    for i, v in enumerate(vals):
        M[i, i] = fmpz(1)
        M[i, n] = fmpz(int(mp.nint(C * v)))
    t0 = time.time()
    R = M.lll(delta=delta, eta=0.51)
    return [[R[r, c] for c in range(n + 1)] for r in range(n)], time.time()-t0

def normalize(row):
    coeffs = [int(c) for c in row[:-1]]
    g = reduce(gcd, [abs(c) for c in coeffs if c])
    coeffs = [c // g for c in coeffs]
    if next(c for c in coeffs if c) < 0:
        coeffs = [-c for c in coeffs]
    return coeffs

mp.mp.dps = 1050
F4 = mp.mpf(open("hex_F4.txt").read())
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
L13, L23, L34 = Li(4, mp.mpf(1)/3), Li(4, mp.mpf(2)/3), Li(4, mp.mpf(3)/4)
l313, l31m = Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2)
z3 = mp.zeta(3); l213 = Li(2, mp.mpf(1)/3); pi2 = mp.pi**2
basis4_names = ["Li4(1/3)", "Li4(2/3)", "Li4(3/4)",
                "Li3(1/3)*log2", "Li3(1/3)*log3", "Li3(-1/2)*log2", "Li3(-1/2)*log3",
                "zeta3*log2", "zeta3*log3",
                "Li2(1/3)^2", "Li2(1/3)*log2^2", "Li2(1/3)*log2*log3", "Li2(1/3)*log3^2",
                "pi2*log2^2", "pi2*log2*log3", "pi2*log3^2", "pi4",
                "log2^4", "log2^3*log3", "log2^2*log3^2", "log2*log3^3", "log3^4"]
basis4 = [L13, L23, L34,
          l313*l2, l313*l3, l31m*l2, l31m*l3,
          z3*l2, z3*l3,
          l213**2, l213*l2**2, l213*l2*l3, l213*l3**2,
          pi2*l2**2, pi2*l2*l3, pi2*l3**2, mp.pi**4,
          l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4]
cand = {}
for dps in (600, 900):
    rows, dt = find_relation_rows([F4] + basis4, dps)
    print(f"dps {dps}: LLL {dt:.1f}s")
    for r in rows[:2]:
        coeffs = normalize(r)
        resid = sum(c*v for c, v in zip(coeffs, [F4] + basis4))
        print(f"   height {max(abs(c) for c in coeffs)}, residual {mp.nstr(abs(resid), 3)}")
    cand[dps] = normalize(rows[0])
print("IDENTICAL:", cand[600] == cand[900])
a = cand[900]
a0 = a[0]
print("\nF4 closed form coefficients (rational):")
for i in range(1, len(a)):
    if a[i]:
        print(f"  {Fr(-a[i], a0)!s:>10}  * {basis4_names[i-1]}")
