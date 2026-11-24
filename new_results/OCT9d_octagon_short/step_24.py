import sys, time, itertools
sys.path.insert(0, ".")
from oct9d_common import *
from syminfra import sh2
import oct8
t0=time.time()
# Full shuffle-ideal check for {2,3,5}: sh(e_a x e_b, e_c) + cubes
cols = []
for a,b,c in itertools.product((2,3,5),(2,3,5),(2,3,5)):
    S = sh2({(a,b): 1}, {(c,): 1})
    cols.append(SS235.coord(S))
for tr in itertools.combinations_with_replacement((2,3,5),3):
    cols.append(col_logprod(tr, SS235))
J = sp.Matrix.hstack(*cols)
print("full shuffle ideal dim ({2,3,5}) =", J.rank(), f"({time.time()-t0:.1f}s)")
# and Li2 symbols (1-w) x w over Gargs span?
L2cols = []
for g in Gargs:
    S = {}
    for x1,c1 in oct8.vec_of_rational(Fr(1)-Fr(g)).items():
        for x2,c2 in oct8.vec_of_rational(Fr(g)).items():
            S[(x1,x2)] = S.get((x1,x2), Fr(0)) + c1*c2
    v = sp.zeros(9,1)
    idx = {p:i for i,p in enumerate(itertools.product((2,3,5),(2,3,5)))}
    for t,c in S.items(): v[idx[t],0] += sp.Rational(c)
    L2cols.append(v)
print("span of Li2-symbols (1-w)xw dim =", sp.Matrix.hstack(*L2cols).rank(), "of 9")
