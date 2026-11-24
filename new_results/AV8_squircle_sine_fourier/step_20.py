import sympy as sp
p, y = sp.symbols('p y', positive=True)
g = -sp.log(1-y**p)/p
H = sp.exp(g)
H3 = sp.simplify(sp.diff(H,p,3).subs(p,2))
# express in terms of s2=y^2, c2=1-y^2, ls=log(y), lc2=log(1-y^2)
H3e = sp.expand(H3*(1-y**2)**sp.Rational(7,2)*64)
H3e = H3e.expand()
print(H3e)
