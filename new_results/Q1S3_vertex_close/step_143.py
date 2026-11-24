import pickle, itertools, time
from fractions import Fraction as Fr
import sympy as sp
B4=('a','p2','q1','q2')
B2=list(itertools.product(B4,repeat=2))
sols=pickle.load(open('sunit_li3.pkl','rb'))
names=pickle.load(open('q3_names.pkl','rb'))['names']
xnew=pickle.load(open('q3_xnew.pkl','rb'))['xnew']
need=sorted(set((v[1],v[2]) for (k,v),e in zip(names,xnew) if k=='dec' and e!=0))
print("needed (v,w) pairs:", need)
def v2_coord(vecz, vec1z):
    v=sp.zeros(16,1)
    for k1,c1 in vec1z.items():
        for k2,c2 in vecz.items():
            v[B2.index((k1,k2))] += sp.Rational(-c1*c2)
    return v
def height(i):
    return sum(abs(c) for c in list(sols[i]['vec_z'].values())+list(sols[i]['vec_1z'].values()))
cols=sorted(range(493), key=height)
M16=sp.Matrix.hstack(*[v2_coord(sols[i]['vec_z'],sols[i]['vec_1z']) for i in cols])
rref,pivots=M16.rref()
Mpiv=sp.Matrix.hstack(*[v2_coord(sols[i]['vec_z'],sols[i]['vec_1z']) for i in [cols[p] for p in pivots]])
def vw_target(V,W):
    v=sp.zeros(16,1)
    v[B2.index((V,W))]=sp.Rational(1,2)
    v[B2.index((W,V))]=sp.Rational(-1,2)
    return v
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
for (V,W) in need:
    if (V,W) in X: continue
    sol=Mpiv.gauss_jordan_solve(vw_target(V,W))[0]
    terms=[(cols[pivots[k]], sp.Rational(sol[k])) for k in range(len(pivots)) if sol[k]!=0]
    X[(V,W)]=terms
    print((V,W), "->", terms)
pickle.dump({'X':X}, open('q3_Xli2.pkl','wb'))
