import sys, json, math
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from flint import fmpz_mat, fmpz
from math import gcd
from functools import reduce
from fractions import Fraction as Fr
from cor5_sym import W4Quotient, li_sym, hpl_sym, lincomb, lneg
from hpl_series import H_series

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
    return [Li(4, mp.mpf(-2)), Li(4, mp.mpf(1)/3), Li(4, mp.mpf(2)/3), Li(4, mp.mpf(1)/2), Li(4, mp.mpf(1)/9),
            l313*l2, l313*l3, l31m*l2, l31m*l3, z3*l2, z3*l3,
            l213**2, l213*l2**2, l213*l2*l3, l213*l3**2, pi2*l213,
            pi2*l2**2, pi2*l2*l3, pi2*l3**2, mp.pi**4,
            l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4]
B = mk_basis()
Q = W4Quotient(2)
Z = (Fr(0),)*4
def vec(j, val=1):
    v = list(Z); v[j] = Fr(val); return tuple(v)
G_c = {0: vec(1), 1: lneg(vec(2)), -1: vec(3)}
letter_map = [(Fr(1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(0)), (Fr(-1), Fr(1))]
sols = json.load(open("sunit23_solutions.json"))
V23 = {Fr(xs): (Fr(fx[1]), Fr(fx[2])) for xs, fx, fy in sols}
OM23 = {Fr(xs): (Fr(fy[1]), Fr(fy[2])) for xs, fx, fy in sols}
li4_args = [Fr(-2), Fr(1,3), Fr(2,3), Fr(1,2), Fr(1,9)]

N = int((1120 + 15) / abs(math.log10(0.5))) + 10
c = mp.mpf(1)/2
for w in [(0,0,-1,1), (-1,1,-1,1), (0,1,-1,1)]:
    t = H_series(w, c, N)
    rows = find_relation_rows([t] + B, 1000)
    a = normalize(rows[0])
    mp.mp.dps = 1120
    resid = sum(cc*v for cc, v in zip(a, [t] + B))
    # numeric projection
    tot = [Fr(0)]*3
    for i, x in enumerate(li4_args):
        cf = Fr(-a[i+1], a[0])
        v = Q.coords(li_sym(4, V23[x], OM23[x]))
        for j in range(3): tot[j] += cf*v[j]
    # pipeline projection
    pipe = [Fr(0)]*3
    for word, c0 in hpl_sym(w, G_c).items():
        nw = tuple(lincomb(l, letter_map) for l in word)
        qc = Q.coords({nw: c0})
        for j in range(3): pipe[j] += qc[j]
    print(f"H{w}: resid {mp.nstr(abs(resid),3)} | numeric proj {tuple(tot)} | pipeline proj {tuple(pipe)}")
