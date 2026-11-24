import sys, time, itertools, pickle, random
sys.path.insert(0, ".")
from oct9d_common import *
from sympy.matrices.normalforms import smith_normal_form
t0=time.time()
# product columns as integer vectors in 27-dim
pcols = []
for g in Gargs:
    for p in (2,3,5): pcols.append(col_prod(sp.Rational(g), p, SS235))
for tr in itertools.combinations_with_replacement((2,3,5),3): pcols.append(col_logprod(tr, SS235))
A = sp.Matrix.hstack(*pcols)
Dglob = sp.ilcm(*[x.q for x in A])
Aint = sp.Matrix([sp.Rational(x)*Dglob for x in A])
snf = smith_normal_form(Aint)
diag = [int(snf[i,i]) for i in range(27) if snf[i,i] != 0]
import math
print("rank:", len(diag), "Smith invariants product:", math.prod(diag))
print("diag:", diag, f"({time.time()-t0:.1f}s)")
