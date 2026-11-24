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
    g = reduce(gcd, [abs(c) for c in coeffs if c])
    coeffs = [c // g for c in coeffs]
    if next((c for c in coeffs if c), 1) < 0: coeffs = [-c for c in coeffs]
    return coeffs

mp.mp.dps = 900
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

def reduce_value(target, dps1=500, dps2=700, work=900):
    cands = {}
    for dps in (dps1, dps2):
        rows = find_relation_rows([target] + mk_basis(), dps)
        cands[dps] = normalize(rows[0])
    mp.mp.dps = work
    coeffs = cands[dps2]
    resid = sum(c*v for c, v in zip(coeffs, [target] + mk_basis()))
    same = cands[dps1] == cands[dps2]
    h = max(abs(c) for c in coeffs)
    return same, h, resid, coeffs

for tag, t in [("Li4(1/4)", Li(4, mp.mpf(1)/4)), ("Li4(-1/2)", Li(4, -mp.mpf(1)/2)),
               ("Li4(8/9)", Li(4, mp.mpf(8)/9)), ("Li4(-1/8)", Li(4, -mp.mpf(1)/8)),
               ("Li4(4/3)", Li(4, mp.mpf(4)/3)), ("Li4(3/2)", Li(4, mp.mpf(3)/2))]:
    same, h, resid, coeffs = reduce_value(t)
    ok = same and abs(resid) < mp.mpf(10)**(-800)
    print(f"{tag}: identical {same}, height {h}, resid {mp.nstr(abs(resid),3)} {'REDUCES' if ok else 'NO'}")
    if ok and tag in ("Li4(1/4)", "Li4(-1/2)"):
        a0 = coeffs[0]
        print("   ", [(B4_NAMES[i-1], str(Fr(-coeffs[i], a0))) for i in range(1, len(coeffs)) if coeffs[i]][:8])
