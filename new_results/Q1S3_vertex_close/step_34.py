import mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
mp.mp.dps=30
gl=GaussLegendre(mp.mp)
nodes = gl.calc_nodes(3, 30)
print(type(nodes), len(nodes), nodes[:2])
