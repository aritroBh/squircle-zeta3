import pickle, sympy as sp, sys
wd = '<project>/research/08_theory/TOW5R_review'
sys.path.insert(0, wd+'/tow5_copy')
d, m = sp.symbols('d m')
a9s = pickle.load(open(wd+'/tow5_copy/a9_closed.pkl','rb'))
a10s = pickle.load(open(wd+'/tow5_copy/a10_closed.pkl','rb'))
a9 = sp.sympify(a9s, locals={'d':d,'m':m,'I':sp.I})
a10 = sp.sympify(a10s, locals={'d':d,'m':m,'I':sp.I})
import tow5lib as T5
A9 = T5.arho(9); A10 = T5.arho(10)
print('a9_closed == arho(9):', sp.cancel(sp.together(a9 - A9)) == 0)
print('a10_closed == arho(10):', sp.cancel(sp.together(a10 - A10)) == 0)
# check the displayed denominator/numerator structure
n9, d9 = sp.fraction(sp.cancel(A9))
n10, d10 = sp.fraction(sp.cancel(A10))
print('a9 den:', sp.factor(d9))
print('a10 den:', sp.factor(d10))
P9 = sp.Poly(sp.expand(n9), d, m)
print('a9 num: deg_d', sp.degree(n9, d), 'terms', len(P9.terms()), 'lc', sp.LC(sp.Poly(n9, d)))
P10 = sp.Poly(sp.expand(n10), d, m)
print('a10 num: deg_d', sp.degree(n10, d), 'terms', len(P10.terms()), 'lc', sp.LC(sp.Poly(n10, d)))
