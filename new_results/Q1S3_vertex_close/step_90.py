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
    Gp.append((height(i), len(Gp), (i,j), coord(g)))
Gp.sort()
cols=[g for _,_,_,g in Gp]+[g for _,g in D]
names=[pr for _,_,pr,_ in Gp]+[d for d,_ in D]
t0=time.time()
M=sp.Matrix.hstack(*cols)
rref,pivots=M.rref()
Mpiv=sp.Matrix.hstack(*[cols[p] for p in pivots])
print("rref+pivots:", len(pivots), " %.0fs"%(time.time()-t0))
S=pickle.load(open('q3_symbols.pkl','rb'))
result={}
for nm,t in [("V3p",coord(S['S3_V3'][0])),("V1A",coord(S['S3_V1'][0])),("V1B",coord(S['S3_V1'][1]))]:
    sol=Mpiv.gauss_jordan_solve(t)[0]
    li3=[(names[p], sp.Rational(sol[k])) for k,p in enumerate(pivots) if p<len(Gp) and sol[k]!=0]
    dec=[(names[p], sp.Rational(sol[k])) for k,p in enumerate(pivots) if p>=len(Gp) and sol[k]!=0]
    result[nm]={'li3':li3,'dec':dec}
    print(nm, ": ReLi3 terms:", len(li3), " decomp terms:", len(dec))
    for pr,c in li3: print("   ReLi3", pr, c)
pickle.dump(result, open('q3_solve.pkl','wb'))
print("total %.0fs"%(time.time()-t0))
