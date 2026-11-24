import pickle, sympy as sp
from o6util import load
mons = load('E3_monomials.pkl')['mons']
t = sp.symbols('t'); r2 = sp.sqrt(2)
v = sp.sqrt(10-7*r2); u = sp.sqrt(2+r2)
# keys classification
from collections import Counter
cls = Counter()
for (a,b,c1,c2,e1,e2), coef in mons:
    kind = 'local' if (c1==0 and c2==0) else ('single1' if c2==0 and c1>0 else ('single2' if c1==0 and c2>0 else 'MIXED'))
    cls[kind]+=1
print(cls)
print("\nMIXED monomials (a,b,c1,c2,e1,e2):")
for key, coef in mons:
    a,b,c1,c2,e1,e2 = key
    if c1>0 and c2>0:
        tpoly = sp.Poly(coef, t)
        print(key, 't-degrees:', [tt[0] for tt in tpoly.terms()])
