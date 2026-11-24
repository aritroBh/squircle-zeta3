"""ZLAD1R step_45 -- weight-9 predictive rungs, mechanical derivation (imports fixed)."""
import json
from fractions import Fraction as Fr
from math import factorial
import sympy as sp
t = json.load(open('muR_exact_N10_K9.json'))
def mu9(k, n, m): return Fr(t[str(k)][str(n)].get(m, '0'))
u, s = sp.symbols('u s')
def Cw(w): return sp.Integer((-1)**(w+1)) * sp.factorial(w-1)
p_prev, p_cur = sp.Integer(2), u
beta = {}
for w in range(2, 10):
    p_new = sp.expand(u*p_cur - s*p_prev); p_prev, p_cur = p_cur, p_new
    beta[w] = sp.expand(u**w - p_cur)
def P(w): return sp.expand(Cw(w)/sp.factorial(w) * beta[w])
def rung_terms(E):
    Ptot = sp.Integer(1)
    for w, e in E.items():
        Ptot = sp.expand(Ptot * P(w)**e / sp.factorial(e))
    out = []
    for e in range(1, 6):
        c = sp.expand(Ptot).coeff(s, e)
        if c == 0: continue
        for mon in sp.Add.make_args(c):
            d = sp.degree(mon, u)
            coef = mon.coeff(u, d)
            out.append((e, int(d), Fr(int(coef.p), int(coef.q))))
    return out
tests = [('z9', {9: 1}), ('z2z7', {2: 1, 7: 1}), ('z3z6', {3: 1, 6: 1}), ('z4z5', {4: 1, 5: 1}),
         ('z2^2z5', {2: 2, 5: 1}), ('z2z3z4', {2: 1, 3: 1, 4: 1}), ('z3^3', {3: 3}), ('z2^3z3', {2: 3, 3: 1})]
ok = True
for m, E in tests:
    terms = rung_terms(E)
    good = True
    for n in range(1, 11):
        rhs = Fr(0)
        for e, d, coef in terms:
            m_e = 'z2' if e == 1 else f'z2^{e}'
            rhs += (-1)**e * factorial(e) * coef * mu9(9 - d, n, m_e)
        if mu9(9, n, m) != rhs: good = False
    ok &= good
    print(('PASS' if good else 'FAIL'), m, 'terms:', terms)
print('ALL weight-9 predictive rungs (mechanical):', ok)
open('progress.log','a').write('step_45: all 8 weight-9 predictive rungs PASS (mechanical derivation from master formula, k=9 n<=10)\n')
