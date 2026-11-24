import sys, time, itertools, pickle
import numpy as np
sys.path.insert(0, ".")
from oct9d_common import *
from sympy.matrices.normalforms import smith_normal_form
t0=time.time()
TW = TW3()
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
pcols = []
for g in Gargs:
    for p in (2,3,5): pcols.append(col_prod(g, p, SS235))
for tr in itertools.combinations_with_replacement((2,3,5),3): pcols.append(col_logprod(tr, SS235))
Pmat235 = sp.Matrix.hstack(*pcols)
cls235, rkP235 = quotient(Pmat235, 27)
tau = {e: cls235(TW[e]) for e in (1,2,3)}
gamQ = {g: cls235(col_Li3(sp.Rational(g), SS235)) for g in Gargs}
gamQ = {g: v for g, v in gamQ.items() if any(v)}

# global denom clear -> SNF
allv = list(gamQ.values())
Dglob = sp.ilcm(*[x.q for v in allv for x in v])
Aint = sp.Matrix.hstack(*[sp.Matrix([sp.Rational(x)*Dglob for x in v]) for v in allv])
print("integer matrix", Aint.rows, "x", Aint.cols, f"({time.time()-t0:.1f}s)")
snf = smith_normal_form(Aint)
diag = [snf[i,i] for i in range(8)]
print("Smith invariants:", diag, f"({time.time()-t0:.1f}s)")
import math
prod = math.prod(int(x) for x in diag)
print("product =", prod, "= index of scaled lattice in Z^8")
print("=> minimal possible det (rational coords) =", sp.Rational(prod, Dglob**8))
pickle.dump({'snf_diag': [int(x) for x in diag], 'Dglob': int(Dglob)}, open("snf235.pkl","wb"))
