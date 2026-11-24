import sys, time, itertools
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()
TW = TW3()

# Product span P for {2,3,5}: all w from Gargs (99) x p in {2,3,5}, plus 10 log^3 symmetrizers
pcols = []
for g in Gargs:
    for p in (2,3,5):
        pcols.append(col_prod(g, p, SS235))
trip3 = list(itertools.combinations_with_replacement((2,3,5),3))
for tr in trip3: pcols.append(col_logprod(tr, SS235))
Pmat = sp.Matrix.hstack(*pcols)
rkP = Pmat.rank()
print(f"{{2,3,5}}: P rank = {rkP} (from {len(pcols)} product columns)  ({time.time()-t0:.1f}s)")

Gmat = sp.Matrix.hstack(*[col_Li3(g, SS235) for g in Gargs])
rkG = Gmat.rank()
rkGP = sp.Matrix.hstack(Gmat, Pmat).rank()
print(f"G rank = {rkG}, G+P rank = {rkGP}, so dim(G∩P) = {rkG+rkP-rkGP}")
# target classes in V/P
for e in (1,2,3):
    r_with = sp.Matrix.hstack(Pmat, sp.Matrix(TW[e])).rank()
    print(f"edge {e}: t_e mod P has dim {r_with - rkP}")
# all three together
T123 = sp.Matrix.hstack(sp.Matrix(TW[1]), sp.Matrix(TW[2]), sp.Matrix(TW[3]))
print("dim span(t1,t2,t3) mod P =", sp.Matrix.hstack(Pmat, T123).rank() - rkP)
