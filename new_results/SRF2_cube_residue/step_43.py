import sympy as sp
a = sp.symbols('a', positive=True)
tt = sp.sqrt(2)*a/sp.sqrt(1+a**4)
integrand_t = sp.simplify(tt**2/sp.sqrt(1-tt**4) * sp.diff(tt, a))
integrand_t = integrand_t.subs(sp.sqrt(a**8 - 2*a**4 + 1), 1-a**4)
print("(iii) simplified:", sp.simplify(integrand_t))
print("equals target:", sp.simplify(integrand_t - 2*sp.sqrt(2)*a**2/(1+a**4)**sp.Rational(3,2))==0)
