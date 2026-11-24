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
ML3=sp.Matrix(27,24,lambda i,j: sp.Rational(LI3v[j][i]))
MD=sp.Matrix(27,len(Dgv),lambda i,j: sp.Rational(Dgv[j][i]))
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
def solve_comp(k, fv):
    tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[k][i]))
    sol=M.gauss_jordan_solve(tgt)[0]
    free=set()
    for z in range(sol.rows): free|=sol[z].free_symbols
    sub={fs:sp.Integer(fv) for fs in free}
    ck=[Fr(sp.Rational(sol[z].subs(sub))) for z in range(24)]
    resid = ML3*sp.Matrix(24,1,lambda i,j: sp.Rational(ck[i])) - tgt
    ok = MD.row_join(resid).rank()==MD.rank()
    return ck, ok
delta={}
for k in range(3):
    ck,ok = solve_comp(k, 0)
    nzc=sum(1 for x in ck if x!=0)
    print(f"component {k}: resid-in-Decomp={ok}, nonzero c_z count={nzc}")
    for z in range(24):
        if ck[z]:
            for key,cf in DL[z].items():
                comp=(cf,Fr(0),Fr(0),1) if k==0 else (Fr(0),cf,Fr(0),1) if k==1 else (Fr(0),Fr(0),cf,1)
                delta[key]=E.fadd(delta.get(key,E.ZERO), comp)
delta={k:v for k,v in delta.items() if not E.feq0(v)}
print("delta terms:", len(delta))
for k,v in sorted(delta.items()): print("  ",k,"=",v[:3],"/",v[3]," ~",E.fev(v,18))
pickle.dump(delta, open('delta_FINAL.pkl','wb'))
open('progress.log','a').write('PG6 s_C7: delta recomputed with resid checks, %d terms\n'%len(delta))
