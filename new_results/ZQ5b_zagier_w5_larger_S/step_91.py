import json, math
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
F = json.load(open('d17_final.json'))
print("T_coef:", F['T_coef'], " z5_coef:", F['z5_coef'])
li5 = F['li5']; prods = F['products']
print("n li5:", len(li5), " n prods:", len(prods))
# sort li5 by |coef| desc for readability? keep insertion order
terms = []
for k,v in li5.items():
    terms.append((k, v))
print("LI5 TERMS (arg : coeff):")
for k,v in terms:
    print(f"  Li5({k}) : {v}")
print("PRODUCTS:")
for k,v in prods.items():
    print(f"  {k} : {v}")
# lcm denominator of the normalized form (divide everything by T_coef=-1419857=-17^5)
from fractions import Fraction as Ff
dens = [Ff(v).denominator for v in list(li5.values())+list(prods.values())] + [Ff(F['z5_coef']).denominator]
import sympy as sp
print("lcm den of integer coeffs:", math.lcm(*dens))
print("17^5 =", 17**5)
