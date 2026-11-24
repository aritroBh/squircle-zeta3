import pickle, itertools, time
from fractions import Fraction as Fr
import sympy as sp
B4=('a','p2','q1','q2')
B2=list(itertools.product(B4,repeat=2))
sols=pickle.load(open('sunit_li3.pkl','rb'))
def v2_coord(vecz, vec1z):
    """S(Li2(z)) = -(1-z) x z as 16-dim vector."""
    v=sp.zeros(16,1)
    for k1,c1 in vec1z.items():
        for k2,c2 in vecz.items():
            v[B2.index((k1,k2))] += sp.Rational(-c1*c2)
    return v
t0=time.time()
M16=sp.Matrix.hstack(*[v2_coord(s['vec_z'],s['vec_1z']) for s in sols])
print("rank[493 Li2 syms]:", M16.rank(), "of 16  (%.0fs)"%(time.time()-t0))
# target: (1/2)(V x W - W x V) for the three pairs
def vw_target(V,W):
    v=sp.zeros(16,1)
    v[B2.index((V,W))]=sp.Rational(1,2)
    v[B2.index((W,V))]=sp.Rational(-1,2)
    return v
X={}
for (V,W) in [('a','p2'),('a','q1'),('a','q2')]:
    t=vw_target(V,W)
    Mx=sp.Matrix.hstack(M16,t)
    print((V,W),"in Li2 span:", Mx.rank()==M16.rank())
    sol=M16.gauss_jordan_solve(t)[0]
    terms=[(i,sp.Rational(sol[i])) for i in range(493) if sol[i]!=0]
    X[(V,W)]=terms
    print("  Li2 terms:", len(terms), " coeffs:", [(i,c) for i,c in terms][:10])
pickle.dump(X, open('q3_Xli2.pkl','wb'))
print("time %.0fs"%(time.time()-t0))
