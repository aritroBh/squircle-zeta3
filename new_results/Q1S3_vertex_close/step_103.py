import pickle, itertools, time
from fractions import Fraction as Fr
import sympy as sp
from q1s_symbol import tadd, t_conj
from q1s_exact import fconj, feq
B4=('a','p2','q1','q2')
BIDX=list(itertools.product(B4,repeat=3))
def coord(T):
    v=sp.zeros(64,1)
    for k,c in T.items(): v[BIDX.index(k)]=sp.Rational(c)
    return v
def deco_gens():
    gens=[]
    for u in B4:
        for v in B4:
            for w in B4:
                t={}
                for perm in [(u,v,w),(v,u,w),(v,w,u)]: t[perm]=t.get(perm,Fr(0))+1
                gens.append(((u,v,w),coord(t)))
    return gens
D=deco_gens()
sols=pickle.load(open('sunit_li3.pkl','rb'))
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
    Gp.append((height(i), (i,j), coord(g)))
Gp.sort(key=lambda x:x[0])
# select independent Re-gens greedily (low height first), then decomp
t0=time.time()
sel=[]; M=sp.zeros(64,0)
for h,pr,g in Gp:
    Mx=sp.Matrix.hstack(M,g)
    if Mx.rank()>M.rank():
        sel.append(('li3',pr)); M=Mx
    if M.rank()>=24: break
print("Re indep selected:", len(sel), "rank:", M.rank())
for d,g in D:
    Mx=sp.Matrix.hstack(M,g)
    if Mx.rank()>M.rank():
        sel.append(('dec',d)); M=Mx
print("total columns:", len(sel), "rank:", M.rank(), " %.0fs"%(time.time()-t0))
pickle.dump({'sel':sel}, open('q3_sel.pkl','wb'))
