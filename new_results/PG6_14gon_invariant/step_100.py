import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import sympy as sp
g=S.init()
letters=[0,1,2]; basis3=list(itertools.product(letters,repeat=3)); bidx={t:i for i,t in enumerate(basis3)}
Dg=[]
for a in letters:
    for b in letters:
        for c in letters:
            v={}
            for t in [(a,b,c),(b,a,c),(b,c,a)]: v[t]=v.get(t,0)+1
            Dg.append(v)
def t2v(T):
    v=[Fr(0)]*27
    for t,c in T.items():
        if t in bidx: v[bidx[t]]+=Fr(c)
    return v
L=pickle.load(open('letters14.pkl','rb')); CRs=[e[0] for e in L['CRset']]
LI3v=[t2v(S.li3_symbol(z)) for z in CRs]; Dgv=[t2v(T) for T in Dg]
# pick a MINIMAL independent subset of the 24 li3 that spans (rank 2 mod Decomp... but need to capture sigma)
# sigma is in the 2-dim span; pick 8 independent li3 to form a clean basis. Find independent columns.
M=sp.Matrix(27,len(LI3v)+len(Dgv),lambda i,j: sp.Rational((LI3v+Dgv)[j][i]))
ML3=sp.Matrix(27,24,lambda i,j: sp.Rational(LI3v[j][i]))
MD=sp.Matrix(27,len(Dgv),lambda i,j: sp.Rational(Dgv[j][i]))
# find a set of li3 indices that are independent mod Decomp
import itertools as it
def rank_with(idxs):
    cols=[LI3v[z] for z in idxs]+Dgv
    return sp.Matrix(27,len(cols),lambda i,j: sp.Rational(cols[j][i])).rank()
# greedy: build independent li3 set
indep=[]
r0=MD.rank()
for z in range(24):
    if rank_with(indep+[z])>rank_with(indep): indep.append(z)
print("independent li3 indices (basis):", indep, " rank:", rank_with(indep)-r0)
# delta_li3
def wedge(u,v):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            if a==b: continue
            key,sg=((a,b),1) if a<b else ((b,a),-1)
            out[key]=out.get(key,Fr(0))+sg*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_li3(z):
    vz=S.vec(z); vomz=S.vec(E.fsub(E.ONE,z)); w=wedge(vz,vomz); out={}
    for ab,cab in w.items():
        for l3,c3 in vz.items(): out[(ab,l3)]=out.get((ab,l3),Fr(0))+cab*c3
    return {k:c for k,c in out.items() if c}
DL=[delta_li3(z) for z in CRs]
sigma=pickle.load(open('sigma14.pkl','rb'))
s_comp=[[Fr(0)]*27 for _ in range(3)]
for key,cFT in sigma.items():
    A,B,C,d=cFT
    if key in bidx: s_comp[0][bidx[key]]+=Fr(A)/d; s_comp[1][bidx[key]]+=Fr(B)/d; s_comp[2][bidx[key]]+=Fr(C)/d
# reduce to MINIMAL basis (unique solution now)
Mmin=sp.Matrix(27,len(indep)+len(Dgv),lambda i,j: sp.Rational(([LI3v[z] for z in indep]+Dgv)[j][i]))
dcomp=[{},{},{}]
for k in range(3):
    tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[k][i]))
    sol=Mmin.gauss_jordan_solve(tgt)[0]
    free=set()
    for z in range(sol.rows): free|=sol[z].free_symbols
    ck=[Fr(sp.Rational(sol[z].subs({fs:0 for fs in free}))) for z in range(len(indep))]
    for zi,z in enumerate(indep):
        if ck[zi]:
            for key,cf in DL[z].items(): dcomp[k][key]=dcomp[k].get(key,Fr(0))+cf*ck[zi]
allk=set(dcomp[0])|set(dcomp[1])|set(dcomp[2])
dmin={}
for key in allk:
    A=dcomp[0].get(key,Fr(0)); B=dcomp[1].get(key,Fr(0)); C=dcomp[2].get(key,Fr(0))
    if A or B or C: dmin[key]=(A,B,C,1)
print("delta from MINIMAL basis:", {k:v[:3] for k,v in sorted(dmin.items())})
dmax=pickle.load(open('delta_FINAL.pkl','rb'))
print("delta from 24-CR basis  :", {k:v[:3] for k,v in sorted(dmax.items())})
print("MATCH:", dmin==dmax)
open('progress.log','a').write('PG6 s_C10: minimal-basis delta matches 24-CR delta: %s\n'%(dmin==dmax))
