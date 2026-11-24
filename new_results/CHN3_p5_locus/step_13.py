setup_code = '''"""common.py: shared setup for CHN3 steps."""
import pickle, sys, os
from fractions import Fraction as F
import sympy as sp

a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
A = [a1,a2,a3,a4,a5]
SA = a1+a2+a3+a4+a5

def dict_to_sym(p, nvars=5):
    vs = A[:nvars]
    return sp.expand(sum(sp.Integer(c.numerator)/c.denominator * sp.prod(vs[i]**e[i] for i in range(nvars)) for e,c in p.items()))

def sym_to_dict(poly, nvars=5):
    P = sp.Poly(sp.expand(poly), *A[:nvars])
    out = {}
    for mon, c in P.terms():
        out[tuple(int(k) for k in mon)] = F(int(c.p), int(c.q))
    return out

def load_forms():
    return pickle.load(open('forms10.pkl','rb'))

def get_C4_Q7(forms=None):
    if forms is None: forms = load_forms()
    P8 = dict_to_sym(forms[8][3])
    C4 = sp.expand(10*P8 / (a2*a3*a4*SA))
    P10 = dict_to_sym(forms[10][6])
    Q7 = sp.expand(224*P10 / (a2*a3*a4*SA))
    return C4, Q7

def peval(p, pt):
    return sum(c*sp.prod(sp.Integer(pt[i])**e[i] for i in range(len(e))) for e,c in p.items())
'''
open('common.py','w').write(setup_code)
print("common.py written")
