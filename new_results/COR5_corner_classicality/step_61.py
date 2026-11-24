import mpmath as mp
from flint import fmpz_mat, fmpz
from fractions import Fraction as Fr
import time

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

mp.mp.dps = 1050
W3 = mp.mpf(open("hex_W3.txt").read())
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
basis3_names = ["zeta3", "Li3(1/3)", "Li3(-1/2)", "pi2*log2", "pi2*log3",
                "Li2(1/3)*log2", "Li2(1/3)*log3", "log2^3", "log2^2*log3", "log2*log3^2", "log3^3"]
basis3 = [mp.zeta(3), Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2), mp.pi**2*l2, mp.pi**2*l3,
          Li(2, mp.mpf(1)/3)*l2, Li(2, mp.mpf(1)/3)*l3, l2**3, l2**2*l3, l2*l3**2, l3**3]

cand = {}
for dps in (600, 900):
    rows = find_relation_rows([W3] + basis3, dps)
    coeffs = rows[0][:-1]
    resid = sum(int(c)*v for c, v in zip(coeffs, [W3] + basis3))
    height = max(abs(int(c)) for c in coeffs)
    print(f"dps {dps}: height {height}, residual {mp.nstr(abs(resid), 3)}")
    cand[dps] = coeffs
print("same relation at both precisions:", cand[600] == cand[900])
a = cand[900]
a0 = int(a[0])
print("exact relation:")
for i in range(1, len(a)):
    if a[i]:
        print(f"  {Fr(-int(a[i]), a0)}  * {basis3_names[i-1]}")
