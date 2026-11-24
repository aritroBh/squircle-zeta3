import sympy as sp
R = sp.Symbol('R')
bh = lambda n: sp.binomial(sp.Rational(1,2), n)
g = (bh(R-1)/bh(R)).rewrite(sp.gamma)
print(repr(g))
print(g.atoms(sp.gamma))
# try simplify with gammasimp
from sympy.simplify.simplify import gammasimp
print(gammasimp(g))
print(sp.simplify(gammasimp(g) - (-2*R/(2*R-3))))
