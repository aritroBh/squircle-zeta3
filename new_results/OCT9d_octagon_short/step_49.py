import sys, time, itertools, pickle
sys.path.insert(0, ".")
from oct9d_common import *
from sympy.matrices.normalforms import smith_normal_form
t0=time.time()
TW = TW3()
pools = pickle.load(open("pools.pkl","rb"))
poolB = pools['poolB']

def colbasis(M):
    rref, pivs = M.rref(); return [M[:,i] for i in pivs]
def quotient(Pmat, dim):
    Bp = colbasis(Pmat); r = len(Bp)
    M = sp.Matrix.hstack(*Bp); E = []
    for i in range(dim):
        ei = sp.zeros(dim,1); ei[i,0]=1
        if (sp.Matrix.hstack(M, *E, ei)).rank() > r+len(E): E.append(ei)
        if r+len(E) == dim: break
    Minv = (sp.Matrix.hstack(M, *E)).inv()
    def cls(v):
        x = Minv*sp.Matrix(v); return x[r:,0]
    return cls, r

# P0: products w in negw7 x p in {2,3,5,7} + 20 log^3
pcols = []
for w in negw7:
    for p in (2,3,5,7): pcols.append(col_prod(sp.Rational(w), p, SS2357))
for tr in itertools.combinations_with_replacement((2,3,5,7),3): pcols.append(col_logprod(tr, SS2357))
Pmat0 = sp.Matrix.hstack(*pcols)
print("P0 rank:", Pmat0.rank(), f"({time.time()-t0:.1f}s)")
cls0, rkP0 = quotient(Pmat0, 64)
print("quotient dim:", 64-rkP0)
tau0 = cls0(TW[0])
# candidate classes: all sols7 (375) with nonzero class
gamB = {}
for g in sols7:
    v = cls0(col_Li3(sp.Rational(g), SS2357))
    if any(v): gamB[g] = v
print("nonzero-class candidates:", len(gamB), f"({time.time()-t0:.1f}s)")
pickle.dump({'tau0': [str(x) for x in tau0]}, open("tau0.pkl","wb"))
# Smith invariants of the class matrix
allv = list(gamB.values())
Dglob = sp.ilcm(*[x.q for v in allv for x in v])
Aint = sp.Matrix.hstack(*[sp.Matrix([sp.Rational(x)*Dglob for x in v]) for v in allv])
print("class matrix", Aint.rows, "x", Aint.cols, f"({time.time()-t0:.1f}s)")
snf = smith_normal_form(Aint)
diag = [int(snf[i,i]) for i in range(20)]
import math
print("Smith invariants:", diag, "product:", math.prod(diag), f"({time.time()-t0:.1f}s)")
pickle.dump({'gamB_keys': list(gamB.keys())}, open("gamB_keys.pkl","wb"))
