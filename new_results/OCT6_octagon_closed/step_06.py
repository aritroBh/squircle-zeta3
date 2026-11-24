import pickle, sympy as sp
from o6util import load
mons = load('E3_monomials.pkl')['mons']
t = sp.symbols('t')
s2 = sp.sqrt(2)
k1 = s2-1; k2 = s2/sp.Integer(2)-1
print("kappa1 =", k1, "=", sp.nsimplify(k1), float(k1))
print("kappa2 =", k2, float(k2))
# Print all 83 monomials with exact coefficients
for i,(key,coeff) in enumerate(mons):
    print(i, key, '|', coeff)
