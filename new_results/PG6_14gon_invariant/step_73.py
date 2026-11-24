import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import sympy as sp
g=S.init()
L=pickle.load(open('letters14.pkl','rb'))
CRs=[e[0] for e in L['CRset']]
LI3=[S.li3_symbol(z) for z in CRs]
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
Dgv=[t2v(T) for T in Dg]; LI3v=[t2v(T) for T in LI3]
Mcols=LI3v+Dgv
M=sp.Matrix(27,len(Mcols),lambda i,j: sp.Rational(Mcols[j][i]))
sigma=pickle.load(open('sigma14.pkl','rb'))
s_comp=[[Fr(0)]*27 for _ in range(3)]
for key,cFT in sigma.items():
    A,B,C,d=cFT
    if key in bidx:
        s_comp[0][bidx[key]]+=Fr(A)/d; s_comp[1][bidx[key]]+=Fr(B)/d; s_comp[2][bidx[key]]+=Fr(C)/d
sol_c=[[Fr(0)]*24 for _ in range(3)]
for k in range(3):
    tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[k][i]))
    sol,params=M.gauss_jordan_solve(tgt, freevar=True)
    sub={p:sp.Integer(0) for p in params}
    for z in range(24):
        v=sol[z]; sol_c[k][z]=Fr(sp.Rational(v.subs(sub)) if v.free_symbols else sp.Rational(v))
Mnull=M.nullspace()
Rli3=[]
for vec in Mnull:
    c=[Fr(sp.Rational(vec[z])) for z in range(24)]
    if any(x!=0 for x in c): Rli3.append(c)
def wedge(u,v):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            if a==b: continue
            key,sg=((a,b),1) if a<b else ((b,a),-1)
            out[key]=out.get(key,Fr(0))+sg*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_li3(z):
    vz=S.vec(z); vomz=S.vec(E.fsub(E.ONE,z))
    w=wedge(vz,vomz); out={}
    for ab,cab in w.items():
        for l3,c3 in vz.items(): out[(ab,l3)]=out.get((ab,l3),Fr(0))+cab*c3
    return {k:c for k,c in out.items() if c}
DL=[delta_li3(z) for z in CRs]
def dadd(*ds):
    o={}
    for d in ds:
        for k,c in d.items(): o[k]=o.get(k,Fr(0))+c
    return {k:c for k,c in o.items() if c}
def dsc(d,q):
    q=Fr(q); return {k:c*q for k,c in d.items() if c*q}
def delta_of_coefs(c):
    acc={}
    for z,cz in enumerate(c):
        if cz: acc=dadd(acc,dsc(DL[z],cz))
    return acc
allzero=all(not delta_of_coefs(r) for r in Rli3)
print("rank[li3|D]=",M.rank()," li3-relations:",len(Rli3)," delta vanishes on relations:",allzero)
delta_sigma={}
for k in range(3):
    dk=delta_of_coefs(sol_c[k])
    for key,c in dk.items():
        comp=(c,Fr(0),Fr(0),1) if k==0 else (Fr(0),c,Fr(0),1) if k==1 else (Fr(0),Fr(0),c,1)
        delta_sigma[key]=E.fadd(delta_sigma.get(key,E.ZERO), comp)
print("delta(sigma) terms:",len(delta_sigma))
for key in sorted(delta_sigma):
    print("  ",key,"=",E.fev(delta_sigma[key],22)," FT",delta_sigma[key])
pickle.dump({'sol_c':sol_c,'delta_sigma':delta_sigma,'Rli3':Rli3,'DL':DL,'CRs':CRs,'allzero':allzero},open('delta14.pkl','wb'))
open('progress.log','a').write('PG6 s_B16: delta(sigma) %d terms, delta on relations vanishes: %s\n'%(len(delta_sigma),allzero))
