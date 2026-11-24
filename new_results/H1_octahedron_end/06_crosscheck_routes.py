"""
H1: cross-check Route A (finite differences of Gamma, uncentred form) against Route B
(cumulants of log-Gamma, centred Dirichlet form) coefficient by coefficient (exact, sympy),
and against the numerical high-order run (05).  Also reduces even zetas to powers of pi.
Usage: python 06_crosscheck_routes.py routeA.json routeB.json hi_coeffs.json n
"""
import sys, json
import sympy as sp
import mpmath as mp
fa, fb, fh, n = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
A = json.load(open(fa)); A = A[n] if n in A else A
B = json.load(open(fb)); B = B[n] if n in B else B
H = json.load(open(fh))
mp.mp.dps = 50
zsub = {sp.Symbol(f'z{2*j}'): sp.zeta(2*j) for j in range(1, 20)}   # even zetas -> pi powers
K = min(len(A['d']), len(B['d'])) - 1
for var in ('eta', 'd'):
    for k in range(K + 1):
        ea, eb = sp.sympify(A[var][k]), sp.sympify(B[var][k])
        diff = sp.expand((ea - eb).subs(zsub))
        num = sp.N(ea.subs({sp.Symbol(f'z{j}'): sp.zeta(j) for j in range(2, 40)}), 45)
        hv = mp.mpf(H[var][k])
        print(f'[{var}^{k}]  A-B (exact, even zetas reduced) = {diff};  value = {sp.N(num, 30)};  |exact - numeric(05)| = {mp.nstr(abs(mp.mpf(str(num)) - hv), 3)}')
