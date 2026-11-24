import mpmath as mp, time
mp.mp.dps = 40
from mpmath.calculus.quadrature import GaussLegendre
q = GaussLegendre(mp.mp)
for d in (1,2,3,4,5,6,7):
    t=time.time()
    r = q.calc_nodes(d, mp.mp.prec)
    print("degree arg", d, "->", len(r), "nodes", f"({time.time()-t:.3f}s)")
