
import mpmath as mp
from mpmath.calculus.quadrature import TanhSinh
mp.mp.dps = 60
TS = TanhSinh(mp.mp)
udeg=5
nodes=[]
for d in range(1, udeg+1): nodes.extend(TS.get_nodes(mp.mpf(0), mp.inf, d, mp.mp.prec))
uh = mp.mpf(2)**(-udeg)
print("num nodes:", len(nodes), " first few (x,w):")
for um,w0 in nodes[:4]: print("   x=%s w=%s"%(mp.nstr(um,10), mp.nstr(w0,10)))
f = lambda x: mp.exp(-x)   # int_0^inf = 1
s1 = sum(w0*uh*f(um) for um,w0 in nodes)
s2 = sum(w0*uh/(4*mp.cosh(um/2)**2)*f(um) for um,w0 in nodes)
print("sum w0*f        =", mp.nstr(s1,20))
print("sum w0/(4cosh^2)*f =", mp.nstr(s2,20))
# test with substitution x = 2 sinh(u/2)? i.e. int f(x) dx = int f(2sinh(u/2)) cosh(u/2) du
s3 = sum(w0*uh*f(2*mp.sinh(um/2))*mp.cosh(um/2) for um,w0 in nodes)
print("sum w0*f(2sinh(u/2))*cosh(u/2) =", mp.nstr(s3,20))
# x = tanh(u/2) maps (0,inf)->(0,1)? no. Maybe get_nodes(0,inf) already for int_0^inf directly?
