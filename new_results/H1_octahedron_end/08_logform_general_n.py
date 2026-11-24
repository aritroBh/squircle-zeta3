"""
H1: 'log form' of the cross-polytope expansion, and the general-n theorem through eta^3.

Claim (PROVED in result.md, Sec. 3): with eta = 1 - 1/p = d/(1+d),
  log[S_n(p)/S_n(1)] = log Rat_n + Lambda_n(eta) + log B_n(eta),
  Rat_n      = (1-eta)^{n-1} / prod_{i=1}^{n-1} (1 - (n-1) eta / i)  ( = prod_{i=1}^{n-2} i/(i-(n-1-i)d) ),
  Lambda_n   = - sum_{j>=2} zeta(j) (n-1)^j (1 - n^{1-j}) eta^j / j          (exact, all orders),
  log B_n    = (n-1)/n zeta(2) eta^2 + 2(n-1)(n+1)/(3n^2) zeta(3) eta^3 + O(eta^4).
Hence the eta^2, eta^3 coefficients of log[S_n/(S_n(1) Rat_n)] are
  -(n-1)(n^2-2n-1)/(2n) zeta(2)   and   -(n^2-1)((n-1)^3-2)/(3n^2) zeta(3).
This script (i) checks those against the exact Route A coefficients for every available n,
(ii) prints log B_n through the available order (the 'bending factor' alone).
"""
import json, sys, glob, re
import sympy as sp

eta = sp.Symbol('eta')
z = {j: sp.Symbol(f'z{j}') for j in range(2, 30)}
files = sorted(glob.glob('routeA_n*_K*.json'))
data = {}
for f in files:
    D = json.load(open(f))
    for n, v in D.items():
        K = len(v['eta']) - 1
        if int(n) not in data or data[int(n)][0] < K:
            data[int(n)] = (K, v['eta'])
even = {z[2 * j]: sp.zeta(2 * j) for j in range(1, 14)}
out = {}
for n in sorted(data):
    K, ce = data[n]
    ser = sum(sp.sympify(c) * eta ** k for k, c in enumerate(ce))
    logS = sp.series(sp.log(ser), eta, 0, K + 1).removeO()
    rat = (1 - eta) ** (n - 1) / sp.prod([1 - sp.Rational(n - 1, i) * eta for i in range(1, n)])
    lograt = sp.series(sp.log(rat), eta, 0, K + 1).removeO()
    Lam = -sum(z[j] * (n - 1) ** j * (1 - sp.Rational(1, n ** (j - 1))) * eta ** j / j for j in range(2, K + 1))
    logB = sp.expand(logS - lograt - Lam)
    pred2 = -sp.Rational((n - 1) * (n * n - 2 * n - 1), 2 * n) * z[2]
    pred3 = -sp.Rational((n * n - 1) * ((n - 1) ** 3 - 2), 3 * n * n) * z[3]
    got2 = sp.expand(logS - lograt).coeff(eta, 2)
    got3 = sp.expand(logS - lograt).coeff(eta, 3)
    print(f'n={n} (K={K}): eta^2 check {sp.simplify((got2 - pred2).subs(even))}, eta^3 check {sp.simplify((got3 - pred3).subs(even))}')
    print(f'   [eta^1] log(S/S1) = {sp.expand(logS).coeff(eta, 1)}   (Rat gives {sp.expand(lograt).coeff(eta, 1)})')
    rows = []
    for k in range(2, K + 1):
        ck = sp.expand(logB.coeff(eta, k))
        ckr = sp.nsimplify(sp.expand(ck.subs(even) / sp.pi ** k)) if all(str(s)[1:].isdigit() and int(str(s)[1:]) % 2 == 0 for s in ck.free_symbols) else None
        rows.append(str(ck))
        print(f'   [eta^{k}] log B_n = {ck}' + (f'   = ({ckr}) pi^{k}' if ckr is not None else ''))
    out[n] = rows
json.dump(out, open('logB_coeffs.json', 'w'), indent=1)
print('saved logB_coeffs.json')
