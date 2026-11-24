
import mpmath as mp, time
from mpmath.calculus.quadrature import TanhSinh
mp.mp.dps = 50
TS = TanhSinh(mp.mp)
t0=time.time()
nodes=[]
for d in range(1, 4):
    nd = TS.get_nodes(mp.mpf(0), mp.inf, d, mp.mp.prec)
    print("deg", d, "nodes:", len(nd), "type", type(nd), "t=%.1f"%(time.time()-t0))
    nodes.extend(nd)
print("total", len(nodes), time.time()-t0)
