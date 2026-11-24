import inspect
from mpmath.calculus import quadrature
print(inspect.getsource(quadrature.QuadratureRule.calc_nodes))
print(inspect.getsource(quadrature.GaussLegendre)[:1500])
