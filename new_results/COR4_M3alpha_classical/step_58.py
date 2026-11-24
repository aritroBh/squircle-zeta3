import sympy as sp
c = sp.Symbol('c')
r = sp.Integer(2)**0 * c**0 * (1 - c)**1 * (1 + c)**(-1)
print("r =", r)
om = sp.expand(1 - r)
print("1-r =", om)
num, den = sp.fraction(om)
print("num, den:", num, den)
print(sp.factor_list(sp.expand(num)))
# check Poly
p = sp.Poly(sp.expand(num), c)
print(p.all_coeffs())
