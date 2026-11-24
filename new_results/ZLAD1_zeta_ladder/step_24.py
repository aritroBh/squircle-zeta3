code = r'''"""step_05.py -- independent verification of the two proof pillars.
P2: (d/da + d/db) E_n(u) = u E_n(u) as formal u-series (symbolic check, n=1..3, order u^6).
P1: [zeta^E]G_n(u) = P^E(d_xi,d_eta) Psi_n(u;0,0) -- checked directly for (n,k) = (2,4), E={2,2},
    by building Psi_n(u;xi,eta) := H-eval of E_n(a+xi,b+eta) explicitly with zeta-free Bell polynomials.
"""
import sympy as sp
from fractions import Fraction as Fr
from math import comb, factorial

u, y, a, b = sp.symbols('u y a b')

def E_series(nv, K):
    l = b - a
    A = -(nv-1)*a + nv*b
    logE = u*A + (1-2*nv)*sp.log(1 + y*(sp.exp(u*l)-1))
    s = sp.series(logE, u, 0, K+1).removeO()
    E = sp.series(sp.exp(s), u, 0, K+1).removeO()
    return sp.expand(E)

print('--- P2: (d_a+d_b)E_n = u E_n ---')
for nv in (1, 2, 3):
    E = E_series(nv, 6)
    D = sp.diff(E, a) + sp.diff(E, b)
    diff = sp.expand(D - u*E)
    print(f'n={nv}: max residual term =', len(sp.Poly(diff, u).terms()) if diff != 0 else 0, '(0 means identity holds)')

print()
print('--- P1 for (n,k)=(2,4), E={2,2}: [z2^2]g_{4,2} = (1/2!) d_xi^2 d_eta^2 Psi_2 [u^4] at 0 ---')
# Build Psi_2(u;xi,eta) = H-eval of E_2(a+xi, b+eta) to order u^4.
# H-eval: y^j a^i b^l -> poch(2,j)/poch(4,j) * DH_{i,l}(2+j, 2),
# DH = Bell in zeta-free cumulants kap^H_{r,s}(p,q):
#   kap^H_{r,s} = (-1)^s psi0^{(r+s-1)}(p)... = harmonic parts:
#   m=r+s=1: H_{p-1} +/- ... ; m>=2: -(-1)^m (m-1)! H^{(m)}_{t-1} pieces (zeta removed).
# We keep harmonic numbers H^{(m)}_t as SYMBOLS to get an exact symbolic value, then evaluate.
xi, eta = sp.symbols('xi eta')
Hsym = {}
def Hh(r, m):  # symbolic harmonic number
    if (r, m) not in Hsym:
        Hsym[(r, m)] = sp.symbols(f'H{r}_{m}')
    return Hsym[(r, m)]

def psiH(m, t):
    # zeta-free part of psi^{(m-1)}(t): m=1 -> H_{t-1}; m>=2 -> -(-1)^m (m-1)! H^{(m)}_{t-1}
    if m == 1:
        return Hh(1, t-1)
    return -(-1)**m * factorial(m-1) * Hh(m, t-1)

def kapH(r, s, p, q):
    m = r+s
    out = psiH(m, p+q)*(-1)
    if s == 0:
        out += psiH(m, p)
    if r == 0:
        out += psiH(m, q)
    return sp.expand(out)

DBH = {}
def DB(i, l, p, q):
    if (i, l) in DBH: return DBH[(i, l)]
    if i == 0 and l == 0:
        res = sp.Integer(1)
    elif i > 0:
        res = sp.expand(sum(comb(i-1, r)*comb(l, s)*kapH(r+1, s, p, q)*DB(i-1-r, l-s, p, q)
                            for r in range(i) for s in range(l+1)))
    else:
        res = sp.expand(sum(comb(l-1, s)*kapH(0, s+1, p, q)*DB(0, l-1-s, p, q) for s in range(l)))
    DBH[(i, l)] = res
    return res

def poch_ratio(nv, j):
    r = Fr(1)
    for t in range(j):
        r *= Fr(nv+t, 2*nv+t)
    return sp.Rational(r.numerator, r.denominator)

nv = 2; K = 4
E = E_series(nv, K)  # in y, a, b
E_sh = sp.expand(E.subs({a: a+xi, b: b+eta}))
# H-evaluate each u-coefficient
Psi = sp.Integer(0)
for (j, i, l), co in sp.Poly(sp.expand(E_sh), y, a, b).terms():
    ex = co.coeff(xi, 0)*1  # keep xi,eta in coeff
    pass
# simpler: iterate over terms of E_sh as poly in (y,a,b,xi,eta)
Psi = sp.Integer(0)
P = sp.Poly(E_sh, y, a, b, xi, eta)
for (j, i, l, ex, et), co in P.terms():
    Psi += co * poch_ratio(nv, j) * DB(i, l, nv+j, nv) * xi**ex * eta**et
Psi = sp.expand(Psi)
# operator (1/2) d_xi^2 d_eta^2, evaluate at 0, take [u^4]
val = sp.Rational(1, 2) * sp.diff(Psi, xi, 2, eta, 2).subs({xi: 0, eta: 0})
g42_z2z2 = sp.expand(val).coeff(u, 4)
# expected: [z2^2] g_{4,2} = -n(2n-1)*[z2^2]mu_{4,2} = -6 * (-2/105) = 4/35
import json
t = json.load(open('mu_exact_N150_K6.json'))
mu42 = Fr(t['4']['2']['z2^2'])
expected = -6*mu42
print('[z2^2]g_{4,2} via operator formula:', g42_z2z2, ' expected:', sp.Rational(expected.numerator, expected.denominator))
print('P1 match:', sp.Rational(expected.numerator, expected.denominator) == g42_z2z2)

# also P1 for E={3} at (n,k)=(2,4): [z3]g_{4,2} = d_xi d_eta (d_xi+d_eta)... operator P_3 = (C_3/3!) square_3
# = (1/3)(D^3 - dx^3 - dy^3)
val3 = sp.Rational(1, 3) * (sp.diff(Psi, xi, 3, eta, 0) + sp.diff(Psi, xi, 0, eta, 3)
                            - sp.diff(Psi, xi, 3) - sp.diff(Psi, eta, 3))
# wait: square_3 = D^3 - dx^3 - dy^3 = (dx+dy)^3 - dx^3 - dy^3 = 3 dx^2 dy + 3 dx dy^2
val3 = sp.Rational(1, 3) * (3*sp.diff(Psi, xi, 2, eta, 1) + 3*sp.diff(Psi, xi, 1, eta, 2)).subs({xi: 0, eta: 0})
g42_z3 = sp.expand(val3).coeff(u, 4)
mu42z3 = Fr(t['4']['2'].get('z3', '0'))
expected3 = -6*mu42z3
print('[z3]g_{4,2} via operator formula:', g42_z3, ' expected:', sp.Rational(expected3.numerator, expected3.denominator),
      ' match:', sp.Rational(expected3.numerator, expected3.denominator) == g42_z3)
'''
open('step_05.py','w').write(code)
import subprocess, os
os.chdir('<project>/research/08_theory/ZLAD1_zeta_ladder')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_05.py'],
                   capture_output=True, text=True, timeout=1800)
print(r.stdout)
print('ERR', r.stderr[-1500:])
