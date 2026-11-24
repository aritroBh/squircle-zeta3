import sys, time, itertools, pickle
import numpy as np
sys.path.insert(0, ".")
from oct9d_common import *
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

# HNF / lattice index of all 96 class vectors
All8 = sp.Matrix.hstack(*[sp.Matrix(v) for v in gamQ.values()])
# Smith normal form via sympy (8 x 96)
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
# clear denominators columnwise to integers (per-column scaling doesn't change lattice span? scaling changes index... 
# use rational HNF: index = product of HNF diagonal of the row-style lattice; easier: compute determinant of the lattice via gcd of 8x8 minors
import math
minors_gcd = None
keys = list(gamQ.keys())
for sel in itertools.combinations(range(len(keys)), 8):
    d = abs(sp.Matrix.hstack(*[sp.Matrix(gamQ[keys[i]]) for i in sel]).det())
    if d == 0: continue
    minors_gcd = d if minors_gcd is None else math.gcd(minors_gcd, d)
    if minors_gcd == 8: break
print("gcd of nonzero 8x8 minors (over sampled) =", minors_gcd, f"({time.time()-t0:.1f}s)")
