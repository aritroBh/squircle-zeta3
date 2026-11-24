import sys, time, itertools
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()
TW = TW3()
# Edge 0: G0 (375 Li3), P0 products
G0 = sp.Matrix.hstack(*[col_Li3(g, SS2357) for g in sols7])
print("G0 rank =", G0.rank(), f"({time.time()-t0:.1f}s)  [max possible 4*10=40]")
pcols = []
for w in negw7:
    for p in (2,3,5,7):
        pcols.append(col_prod(w, p, SS2357))
trip4 = list(itertools.combinations_with_replacement((2,3,5,7),3))
for tr in trip4: pcols.append(col_logprod(tr, SS2357))
P0 = sp.Matrix.hstack(*pcols)
rkP0 = P0.rank()
print("P0 rank =", rkP0, f"({time.time()-t0:.1f}s)")
print("G0+P0 rank =", sp.Matrix.hstack(G0,P0).rank(), " dim G0∩P0 =", G0.rank()+rkP0-sp.Matrix.hstack(G0,P0).rank())
t0v = sp.Matrix(TW[0])
print("t_0 mod P0 dim =", sp.Matrix.hstack(P0, t0v).rank() - rkP0)
print("t_0 in G0 alone?", G0.rank() == sp.Matrix.hstack(G0, t0v).rank())
