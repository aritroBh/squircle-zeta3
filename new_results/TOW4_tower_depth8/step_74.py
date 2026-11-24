import sys, pickle
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import tow4lib as T4
d, m = T4.d, T4.m
a7 = T4.arho(7); a8 = T4.arho(8)
for nm, a in [('a_7',a7),('a_8',a8)]:
    num, den = sp.fraction(a)
    print(nm, "denominator:", sp.factor(den))
    print(nm, "numerator expanded:")
    print(sp.expand(num))
    print()
