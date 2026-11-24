"""
H1, n = 3: exact 'log form'   log[ (2-p) S(p) / (4 sqrt3) ] = sum_k L_k eta^k,  eta = 1 - 1/p,
from the exact Route A eta-coefficients (routeA_n3_K12.json).  (1/(2-p) = 1/(1-d) is the rational factor Rat_3.)
"""
import json
import sympy as sp
eta = sp.Symbol('eta')
z = {j: sp.Symbol(f'z{j}') for j in range(2, 30)}
even = {z[2 * j]: sp.zeta(2 * j) for j in range(1, 14)}
ce = json.load(open('routeA_n3_K12.json'))['3']['eta']
K = len(ce) - 1
ser = sum(sp.sympify(c) * eta ** k for k, c in enumerate(ce))
# (2-p) = (1-2eta)/(1-eta) * p ... use Rat_3 = (1-eta)/(1-2eta) exactly:  S/S(1) = Rat_3 * exp(sum L_k eta^k)
L = sp.series(sp.log(ser * (1 - 2 * eta) / (1 - eta)), eta, 0, K + 1).removeO()
for k in range(K + 1):
    ck = sp.expand(sp.expand(L).coeff(eta, k).subs(even))
    print(f'L_{k} = {ck}   = {sp.N(ck, 25)}')
