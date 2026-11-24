import sys, json
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from flint import fmpz_mat, fmpz
from math import gcd
from functools import reduce
from fractions import Fraction as Fr
from cor5_sym import W4Quotient, li_sym
from hpl_series import Li31

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

mp.mp.dps = 800
l2, l3 = mp.log(2), mp.log(3)
def Li(s, x): return mp.polylog(s, x)
l313, l31m = Li(3, mp.mpf(1)/3), Li(3, -mp.mpf(1)/2)
z3 = mp.zeta(3); l213 = Li(2, mp.mpf(1)/3); pi2 = mp.pi**2
Q = W4Quotient(2)
sols = json.load(open("sunit23_solutions.json"))
V23 = {Fr(xs): (Fr(fx[1]), Fr(fx[2])) for xs, fx, fy in sols}
OM23 = {Fr(xs): (Fr(fy[1]), Fr(fy[2])) for xs, fx, fy in sols}
def mk_basis(carriers):
    return [Li(4, mp.mpf(x.numerator)/x.denominator) for x in carriers] + \
        [l313*l2, l313*l3, l31m*l2, l31m*l3, z3*l2, z3*l3,
         l213**2, l213*l2**2, l213*l2*l3, l213*l3**2, pi2*l213,
         pi2*l2**2, pi2*l2*l3, pi2*l3**2, mp.pi**4,
         l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4]
Hv = -Li31(mp.mpf(-1)/2, mp.mpf(-1))
for carriers in [[Fr(-2), Fr(3,4), Fr(1,2), Fr(1,9)], [Fr(-2), Fr(1,3), Fr(2,3), Fr(1,2), Fr(1,9)]]:
    B = mk_basis(carriers)
    cands = {}
    for dps in (500, 700):
        rows = find_relation_rows([Hv] + B, dps)
        cands[dps] = normalize(rows[0])
    mp.mp.dps = 800
    resid = sum(cc*v for cc, v in zip(cands[700], [Hv] + B))
    same = cands[500] == cands[700]
    tot = [Fr(0)]*3
    a = cands[700]
    for i, x in enumerate(carriers):
        cf = Fr(-a[i+1], a[0])
        v = Q.coords(li_sym(4, V23[x], OM23[x]))
        for j in range(3): tot[j] += cf*v[j]
    print(f"carriers {carriers}: identical {same}, resid {mp.nstr(abs(resid),3)}, class {tuple(tot)}")
