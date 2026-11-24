
import mpmath as mp, time
from mpmath.calculus.quadrature import TanhSinh
mp.mp.dps = 50
TS = TanhSinh(mp.mp)
nodes=[]
for d in range(1, 4): nodes.extend(TS.get_nodes(mp.mpf(0), mp.inf, d, mp.mp.prec))
uh = mp.mpf(2)**(-3)
print("max node x:", mp.nstr(max(um for um,w0 in nodes),6), " min:", mp.nstr(min(um for um,w0 in nodes),6))
f = lambda x: mp.exp(-x)
t0=time.time()
s1 = mp.fsum([w0*uh*f(um) for um,w0 in nodes])
print("sum w0*f        =", mp.nstr(s1,20), "t=%.1f"%(time.time()-t0))
t0=time.time()
s2 = mp.fsum([w0*uh/(4*mp.cosh(um/2)**2)*f(um) for um,w0 in nodes])
print("sum w0/(4cosh^2)*f =", mp.nstr(s2,20), "t=%.1f"%(time.time()-t0))
