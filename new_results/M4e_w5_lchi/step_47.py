import json
from fractions import Fraction as F
M4E = '<project>/research/08_theory/M4e_w5_lchi'
d = json.load(open(M4E+'/lchi_testB_1200.json'))
names = d['names']; ca=[F(c) for c in d['ca']]; cb=[F(c) for c in d['cb']]
print("c0:", d['c0'], " L coef:", d['aL'], "+", d['bL'], "sqrt2")
print("\nLi5 coefficients (rat + sqrt2*sq):")
for i,nm in enumerate(names):
    if nm.startswith('Li5'):
        print(f"  {nm:28s} {ca[i]}  {cb[i]}")
print("\nmax |num| among Li5:", max(max(abs(ca[i].numerator),abs(cb[i].numerator)) for i,nm in enumerate(names) if nm.startswith('Li5')))
print("num nonzero product terms:", sum(1 for i,nm in enumerate(names) if not nm.startswith('Li5') and (ca[i] or cb[i])))
