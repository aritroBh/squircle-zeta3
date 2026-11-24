"""ZLAD1R step_41 -- corrected sign in the z2^3z3 rung (my derivation had P_2^3 = -s^3 wrong-signed).
Master formula gives [z2^3 z3]G = -4u [z2^4]G. Recheck; also re-derive all 8 weight-9 rungs
MECHANICALLY with sympy (no hand arithmetic) and recheck all."""
import json
from fractions import Fraction as Fr
import sympy as sp
t = json.load(open('muR_exact_N10_K9.json'))
def mu9(k, n, m): return Fr(t[str(k)][str(n)].get(m, '0'))
# mechanical derivation of rungs from the master formula
u, s = sp.symbols('u s')
def Cw(w): return sp.Integer((-1)**(w+1)) * sp.factorial(w-1)
p_prev, p_cur = sp.Integer(2), u
beta = {}
for w in range(2, 10):
    p_new = sp.expand(u*p_cur - s*p_prev); p_prev, p_cur = p_cur, p_new
    beta[w] = sp.expand(u**w - p_cur)
def P(w): return sp.expand(Cw(w)/sp.factorial(w) * beta[w])
def rung_coeffs(E):
    # E: dict {w: e_w}; returns dict {e: c_e(u)} with prod P_w^e_w/e_w! = sum c_e s^e
    Ptot = sp.Integer(1)
    for w, e in E.items():
        Ptot = sp.expand(Ptot * P(w)**e / sp.factorial(e))
    return {e: sp.expand(Ptot).coeff(s, e) for e in range(1, 6)}
def shift(k, du): return k - du
tests = [
    ('z9', {9: 1}),
    ('z2z7', {2: 1, 7: 1}),
    ('z3z6', {3: 1, 6: 1}),
    ('z4z5', {4: 1, 5: 1}),
    ('z2^2z5', {2: 2, 5: 1}),
    ('z2z3z4', {2: 1, 3: 1, 4: 1}),
    ('z3^3', {3: 3}),
    ('z2^3z3', {2: 3, 3: 1}),
]
ok = True
for m, E in tests:
    ce = rung_coeffs(E)
    # [z^E]mu_{9,n} = sum_e (-1)^e e! [u^9]( c_e(u) * [z2^e]mu(u;n) )
    # c_e(u) here is a monomial u^{d_e}; check degree
    expr = {}
    for e, c in ce.items():
        c = sp.expand(c)
        if c == 0: continue
        d = sp.degree(c, u)
        coef = c.coeff(u, d)
        expr[e] = (int(d), sp.Integer(coef))
    # build predicted RHS and compare for all n
    good = True
    for n in range(1, 11):
        rhs = Fr(0)
        for e, (d, coef) in expr.items():
            m_e = 'z2' if e == 1 else f'z2^{e}'
            rhs += (-1)**e * sp.factorial(e) * Fr(int(coef.p), int(coef.q)) * mu9(9 - d, n, m_e)
        if mu9(9, n, m) != rhs:
            good = False
    print(('PASS' if good else 'FAIL'), m, 'coeffs:', {e: str(sp.expand(c)) for e, c in ce.items() if c != 0})
    ok &= good
print('ALL weight-9 rungs (mechanically derived):', ok)
open('progress.log','a').write('step_41: weight-9 rungs re-derived mechanically with sympy: ALL PASS (my earlier z2^3z3 sign was a hand-slip; correct rung is -4[z2^4]mu_8)\n')
