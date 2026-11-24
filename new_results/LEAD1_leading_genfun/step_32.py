import json, sympy as sp
from fractions import Fraction as Fr
import math
lam, pi = sp.symbols('lam pi')
forms = {int(w): sp.sympify(v, locals={'lam':lam,'pi':pi,'z3':sp.Symbol('z3'),'z5':sp.Symbol('z5'),'z7':sp.Symbol('z7'),'z9':sp.Symbol('z9')}) for w,v in json.load(open('cw_closed_forms.json')).items()}
lll = json.load(open('lll_relations_dps400.json'))
def monoms(w):
    out = []
    def rec(rem, factors, start):
        if rem == 0: out.append(list(factors)); return
        for g in range(start, rem+1):
            if g in (1,2) or g%2==1: rec(rem-g, factors+[g], g)
    rec(w, [], 1)
    return out
def sym_of_m(m):
    s = sp.Integer(1)
    for g in m: s *= lam if g==1 else (pi**2 if g==2 else sp.Symbol(f'z{g}'))
    return s
ok_all = True
for w in [5,6,7,8,9]:
    basis = monoms(w)
    f = sp.expand(forms[w])
    coeffs = [sp.expand(f).coeff(sym_of_m(m)) for m in basis]  # exact rationals
    den = 1
    for c in coeffs: den = sp.ilcm(den, sp.denom(c))
    ivec = [int(sp.nsimplify(c*den)) for c in coeffs] + [-int(den)]
    g = 0
    for x in ivec: g = math.gcd(g, abs(x))
    ivec = [x//g for x in ivec]
    if w in lll:
        lrel = lll[str(w)]['rel']
        same = ivec == lrel or ivec == [-x for x in lrel]
        ok_all &= same
        print(f"w={w}: LLL relation == law closed form (up to sign): {same}")
    else:
        print(f"w={w}: (control, no LLL) exact vector {ivec}")
print("ALL MATCH:", ok_all)
