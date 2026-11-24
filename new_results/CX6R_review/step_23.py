import mpmath as mp, time
mp.mp.dps = 40
from mpmath.calculus.quadrature import GaussLegendre
t=time.time()
q = GaussLegendre(mp.mp)
x, w = q.calc_nodes(8, mp.mp.prec)
print("deg 8 took", time.time()-t, "s; nodes:", mp.nstr(x[0],8), mp.nstr(w[0],8))
t=time.time()
x, w = q.calc_nodes(48, mp.mp.prec)
print("deg 48 took", time.time()-t)
t=time.time()
x, w = q.calc_nodes(48, 500)
print("deg 48 @500 bits took", time.time()-t)
