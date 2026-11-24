import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *

sols = pickle.load(open('solsF_B8.pkl','rb'))
zs = sorted(sols.keys(), key=pfloat)

# ---- weight-2 classes: w(z) = (1-z) ^ z in Lambda^2 V (basis: 2^3, 2^e, 3^e) ----
def wedge2(z):
    a = vecfun(psub(ONE,z)); b = vecfun(z)
    # a ^ b over basis ('2','3','e')
    out = {}
    for i,u in enumerate(BASIS):
        for j,v in enumerate(BASIS):
            if j<=i: continue
            c = a.get(u,Fr(0))*b.get(v,Fr(0)) - a.get(v,Fr(0))*b.get(u,Fr(0))
            if c: out[(u,v)] = c
    return out
W = []
for z in zs:
    w = wedge2(z)
    W.append(sp.Matrix([[w.get(('2','3'),0)],[w.get(('2','e'),0)],[w.get(('3','e'),0)]]))
# rank of all Li2 classes
M = sp.Matrix.hstack(*W)
print("rank of Li2-classes over all 129 z:", M.rank(), "(expect 3 = dim Lambda^2 V)")
# pick 3 nice reps
def niceness(i):
    z = zs[i]; _,v2,v3,k = factor_sunit(z)
    fz = pfloat(z)
    pref = 0 if 0 < fz < 1 else (1 if fz < 0 else 2)
    return (pref, abs(v2)+abs(v3)+abs(k))
order = sorted(range(len(zs)), key=niceness)
b2 = []
BM = sp.zeros(3,0)
for i in order:
    M2 = sp.Matrix.hstack(BM, W[i]) if BM.cols else W[i]
    if BM.cols==0 or M2.rank()>BM.rank():
        b2.append(i); BM = M2
    if len(b2)==3: break
print("Li2-class reps:")
for i in b2:
    print("   z=%-14s ~%+.8f fac %s"%(str(zs[i]), pfloat(zs[i]), factor_sunit(zs[i])))
pickle.dump({'b2_idx':b2, 'b2_zs':[zs[i] for i in b2]}, open('li2_basis.pkl','wb'))

# check: Li2(x_1), Li2(x_2) classes spanned?
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
for nm,z in [('x1',Xs[1]),('x2',Xs[2]),('x5',Xs[5]),('x6',Xs[6])]:
    w = wedge2(z)
    v = sp.Matrix([[w.get(('2','3'),0)],[w.get(('2','e'),0)],[w.get(('3','e'),0)]])
    print("   class(%s)="%nm, list(v), "in span of reps:", sp.Matrix.hstack(BM,v).rank()==3)
