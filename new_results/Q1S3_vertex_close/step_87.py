import pickle, itertools
from fractions import Fraction as Fr
import sympy as sp
from q1s_symbol import tadd, tscale, t_conj
import time

B4=('a','p2','q1','q2')
BIDX=list(itertools.product(B4,repeat=3))
def coord(T):
    v=sp.zeros(64,1)
    for k,c in T.items(): v[BIDX.index(k)]=sp.Rational(c)
    return v
def deco_gens():
    """u shuffle (v x w) for u,v,w basis letters (ordered)."""
    gens=[]
    for u in B4:
        for v in B4:
            for w in B4:
                t={}
                for perm in [(u,v,w),(v,u,w),(v,w,u)]: t[perm]=t.get(perm,Fr(0))+1
                gens.append(((u,v,w), coord(t)))
    return gens
t0=time.time()
D=deco_gens()
sols=pickle.load(open('sunit_li3.pkl','rb'))
from q1s_exact import fconj, feq
pairs=[]; used=set()
for i,s in enumerate(sols):
    if i in used: continue
    j=next((jj for jj,s2 in enumerate(sols) if jj!=i and feq(s2['z'],fconj(s['z']))), None)
    used.add(i)
    if j is not None: used.add(j)
    pairs.append((i,j))
def height(i):
    return sum(abs(c) for c in list(sols[i]['vec_z'].values())+list(sols[i]['vec_1z'].values()))
Gp=[]
for i,j in pairs:
    g=tadd(sols[i]['S'], t_conj(sols[i]['S']))
    Gp.append((height(i), i, j, coord(g)))
Gp.sort(key=lambda x:x[0])
print("pairs:", len(Gp), " decomp gens:", len(D))
M=sp.Matrix.hstack(*[g for _,_,_,g in Gp]+[g for _,g in D])
print("rank[ReLi3 | Decomp] =", M.rank(), "of 64   (%.1fs)"%(time.time()-t0))
rrefM, pivots = M.rref()
print("rref done, pivots:", len(pivots), " %.1fs"%(time.time()-t0))
pickle.dump({'pivots':pivots}, open('q3_span_tmp.pkl','wb'))
