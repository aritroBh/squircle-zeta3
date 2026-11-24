import pickle, itertools, time
from fractions import Fraction as Fr
import sympy as sp
from q1s_symbol import tadd, t_conj
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
D={d:g for d,g in deco_gens()}
sols=pickle.load(open('sunit_li3.pkl','rb'))
from q1s_exact import fconj, feq
pairs=[]; used=set()
for i,s in enumerate(sols):
    if i in used: continue
    j=next((jj for jj,s2 in enumerate(sols) if jj!=i and feq(s2['z'],fconj(s['z']))), None)
    used.add(i)
    if j is not None: used.add(j)
    pairs.append((i,j))
pairidx={p:k for k,p in enumerate(pairs)}
def height(i):
    return sum(abs(c) for c in list(sols[i]['vec_z'].values())+list(sols[i]['vec_1z'].values()))
Gp={}
for i,j in pairs:
    g=tadd(sols[i]['S'], t_conj(sols[i]['S']))
    Gp[(i,j)]=coord(g)
sel=pickle.load(open('q3_sel.pkl','rb'))['sel']
cols=[Gp[pr] if k=='li3' else D[d] for k,(d) in [(k,v) for k,v in sel] for pr in [d]]
# fix: build cols properly
cols=[]
for kind,val in sel:
    cols.append(Gp[val] if kind=='li3' else D[val])
M=sp.Matrix.hstack(*cols)
S=pickle.load(open('q3_symbols.pkl','rb'))
out={}
t0=time.time()
for nm,t in [("V3p",coord(S['S3_V3'][0])),("V1A",coord(S['S3_V1'][0])),("V1B",coord(S['S3_V1'][1]))]:
    x=M.gauss_jordan_solve(t)[0]
    terms=[(sel[i], sp.Rational(x[i])) for i in range(len(sel)) if x[i]!=0]
    out[nm]=terms
    li3=[(v,c) for (k,v),c in terms if k=='li3']
    dec=[(v,c) for (k,v),c in terms if k=='dec']
    mx=max(abs(c) for _,c in terms)
    print(nm,"terms:",len(terms)," li3:",len(li3)," dec:",len(dec)," max|coeff|:",mx)
pickle.dump(out, open('q3_solve_min.pkl','wb'))
print("time %.0fs"%(time.time()-t0))
for nm in out:
    print(nm, out[nm])
