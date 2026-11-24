import sys
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from flint import fmpz_mat, fmpz
from math import gcd
from functools import reduce
from fractions import Fraction as Fr

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
F4 = mp.mpf(open("hex_F4.txt").read())
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
l313, l31m = Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2)
z3 = mp.zeta(3); l213 = Li(2, mp.mpf(1)/3); pi2 = mp.pi**2
B4_NAMES = ["Li4(1/3)", "Li4(2/3)", "Li4(3/4)", "Li4(1/2)", "Li4(1/9)",
            "Li3(1/3)*log2", "Li3(1/3)*log3", "Li3(-1/2)*log2", "Li3(-1/2)*log3",
            "zeta3*log2", "zeta3*log3",
            "Li2(1/3)^2", "Li2(1/3)*log2^2", "Li2(1/3)*log2*log3", "Li2(1/3)*log3^2", "pi2*Li2(1/3)",
            "pi2*log2^2", "pi2*log2*log3", "pi2*log3^2", "pi4",
            "log2^4", "log2^3*log3", "log2^2*log3^2", "log2*log3^3", "log3^4"]
def mk_basis():
    return [Li(4, mp.mpf(1)/3), Li(4, mp.mpf(2)/3), Li(4, mp.mpf(3)/4), Li(4, mp.mpf(1)/2), Li(4, mp.mpf(1)/9),
            l313*l2, l313*l3, l31m*l2, l31m*l3, z3*l2, z3*l3,
            l213**2, l213*l2**2, l213*l2*l3, l213*l3**2, pi2*l213,
            pi2*l2**2, pi2*l2*l3, pi2*l3**2, mp.pi**4,
            l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4]
B = mk_basis()
cands = {}
for dps in (700, 1000):
    rows = find_relation_rows([F4] + B, dps)
    cands[dps] = normalize(rows[0])
    mp.mp.dps = 1120
    resid = sum(c*v for c, v in zip(cands[dps], [F4] + B))
    print(f"F4 dps {dps}: height {max(abs(c) for c in cands[dps])}, resid@1120 {mp.nstr(abs(resid),3)}")
same = cands[700] == cands[1000]
print("identical:", same)
if same:
    a = cands[1000]; a0 = a[0]
    print("\nF4 =")
    for i in range(1, len(a)):
        if a[i]:
            print(f"  {Fr(-a[i], a0)!s:>10} * {B4_NAMES[i-1]}")
    import json
    json.dump({"a0": a0, "coeffs": a[1:]}, open("F4_relation.json", "w"))
