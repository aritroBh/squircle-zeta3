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
M=sp.Matrix(27,len(LI3v)+len(Dgv),lambda i,j: sp.Rational((LI3v+Dgv)[j][i]))
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
# compute delta per component (plain rational accumulation)
dcomp=[{},{},{}]
for k in range(3):
    tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[k][i]))
    sol=M.gauss_jordan_solve(tgt)[0]
    free=set()
    for z in range(sol.rows): free|=sol[z].free_symbols
    sub={fs:sp.Integer(0) for fs in free}
    ck=[Fr(sp.Rational(sol[z].subs(sub))) for z in range(24)]
    for z in range(24):
        if ck[z]:
            for key,cf in DL[z].items():
                dcomp[k][key]=dcomp[k].get(key,Fr(0))+cf*ck[z]
# assemble FT delta: key -> (dcomp0, dcomp1, dcomp2)
allk=set(dcomp[0])|set(dcomp[1])|set(dcomp[2])
delta={}
for key in allk:
    A=dcomp[0].get(key,Fr(0)); B=dcomp[1].get(key,Fr(0)); C=dcomp[2].get(key,Fr(0))
    if A or B or C: delta[key]=(A,B,C,1)
print("FINAL delta(N3) of the 14-gon (2,1) invariant:")
for k,v in sorted(delta.items()):
    Ln={1:'a',2:'u',0:'2'}
    print(f"  ({Ln[k[0][0]]}^{Ln[k[0][1]]}) x {Ln[k[1]]} : ({v[0]}) + ({v[1]}) a + ({v[2]}) a^2   ~ {E.fev(v,20)}")
print("NONZERO:", len(delta)>0)
pickle.dump(delta, open('delta_FINAL.pkl','wb'))
open('progress.log','a').write('PG6 s_C8: FINAL delta(N3) = %d terms, nonzero=%s\n'%(len(delta),len(delta)>0))
