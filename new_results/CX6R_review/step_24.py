import mpmath as mp, time
mp.mp.dps = 40
from mpmath.calculus.quadrature import GaussLegendre
q = GaussLegendre(mp.mp)
r = q.calc_nodes(8, mp.mp.prec)
print(type(r), len(r))
print(r[:4])
