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
    Gp.append((height(i),(i,j),coord(g)))
Gp.sort(key=lambda x:x[0])
cols=[g for _,_,g in Gp]+[g for _,g in D]
names=[('li3',pr) for _,pr,_ in Gp]+[('dec',d) for _,d in D]
t0=time.time()
M=sp.Matrix.hstack(*cols)
print("M:", M.rows, M.cols)
K=M.nullspace()
print("kernel dim:", len(K), " %.0fs"%(time.time()-t0))
pickle.dump({'names':names}, open('q3_names.pkl','wb'))
# convert kernel to integer matrix
import sympy
Kint=[]
for kv in K:
    den=1
    for e in kv: den=sp.ilcm(den, sp.denom(e))
    Kint.append([int(sp.nsimplify(e*den)) for e in kv])
pickle.dump({'Kint':Kint}, open('q3_kernel.pkl','wb'))
print("kernel integral, max entry:", max(max(abs(e) for e in kv) for kv in Kint), " %.0fs"%(time.time()-t0))
