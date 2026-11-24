import sympy as sp
R = sp.Symbol('R')
bh = lambda n: sp.binomial(sp.Rational(1,2), n)
ratio = bh(R-1)/bh(R)
g = ratio.rewrite(sp.gamma)
print(sp.simplify(g))
print(sp.simplify(g - (-2*R/(2*R-3))))
