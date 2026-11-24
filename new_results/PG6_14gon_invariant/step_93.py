import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import sympy as sp
g=S.init()
# rebuild both sigmas
def lam_dod3(x,y):
    L1=S.vec(E.fsub(E.fmul(x,y),E.ONE)); L2=S.vec(E.fmul(x,y)); L3=S.vec(E.fsub(y,E.ONE))
    L4=S.vec(E.fdiv(E.fmul(x,E.fsub(y,E.ONE)),E.fsub(x,E.ONE))); L5=S.vec(E.fsub(x,E.ONE)); L6=S.vec(x)
    return S.sadd(S.tensor(L1,L2,L3),S.tensor(L1,L4,L6),S.tensor(L3,L5,L6))
def lam_mine(x,y):
    Sdod=lam_dod3(x,y); out={}
    for (a,b,c),cf in Sdod.items():
        out[(c,b,a)]=out.get((c,b,a),Fr(0))-cf
        if out[(c,b,a)]==0: del out[(c,b,a)]
    return out
S14=pickle.load(open('S14.pkl','rb'))
def is_w3(nm): return nm.startswith('L3x') or nm.startswith('Ax') or nm.startswith('Lam') or nm.startswith('T3p') or nm.startswith('T3_')
def atom_sym_nonlam(name):
    if name.startswith('L3x'): return S.sym_L3x(int(name[3:]))
    if name.startswith('Ax'): return S.sym_Ax(int(name[2:]))
    if name.startswith('T3p'): return S.sym_T3p(int(name[4]),int(name[6]))
    if name.startswith('T3_'): return S.sym_T3a(int(name[3]),int(name[5]),int(name[6]))
def build_sigma(lamfun):
    sigma={}
    for mon,c in S14.items():
        syms=[s for s in dict(mon) if s not in ('l2','la','lu','z2','z3')]
        others=[s for s in dict(mon) if s in ('l2','la','lu','z2','z3')]
        if len(syms)==1 and not others and is_w3(syms[0]):
            nm=syms[0]
            if nm.startswith('Lam'): sym=lamfun(E.X[int(nm[3])], E.fdiv(E.X[int(nm[5])],E.X[int(nm[3])])); scalar=E.ONE
            else: sym,scalar=atom_sym_nonlam(nm)
            tot=E.fmul(c,scalar)
            for key,cf in sym.items(): sigma[key]=E.fadd(sigma.get(key,E.ZERO), E.fs(tot,cf))
    return sigma
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
def delta_full(lamfun):
    sigma=build_sigma(lamfun)
    s_comp=[[Fr(0)]*27 for _ in range(3)]
    for key,cFT in sigma.items():
        A,B,C,d=cFT
        if key in bidx: s_comp[0][bidx[key]]+=Fr(A)/d; s_comp[1][bidx[key]]+=Fr(B)/d; s_comp[2][bidx[key]]+=Fr(C)/d
    dsig={}
    for k in range(3):
        tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[k][i]))
        sol=M.gauss_jordan_solve(tgt)[0]
        free=set()
        for z in range(sol.rows): free|=sol[z].free_symbols
        sub={fs:0 for fs in free}
        ck=[Fr(sp.Rational(sol[z].subs(sub))) for z in range(24)]
        for z in range(24):
            if ck[z]:
                for key,cf in DL[z].items():
                    dsig[key]=E.fadd(dsig.get(key,E.ZERO), ((cf,Fr(0),Fr(0),1) if k==0 else (Fr(0),cf,Fr(0),1) if k==1 else (Fr(0),Fr(0),cf,1)))
    return {k:v for k,v in dsig.items() if not E.feq0(v)}
d1=delta_full(lam_mine); d2=delta_full(lam_dod3)
print("delta (Lam -tau):")
for k,v in sorted(d1.items()): print("  ",k,"=",v[:3],"/",v[3],"  eval",E.fev(v,15))
print("delta (Lam dod3-unconv):")
for k,v in sorted(d2.items()): print("  ",k,"=",v[:3],"/",v[3],"  eval",E.fev(v,15))
print("identical:", all(E.feq0(E.fsub(d1.get(k,E.ZERO),d2.get(k,E.ZERO))) for k in set(d1)|set(d2)))
pickle.dump({'d1':d1,'d2':d2}, open('delta_final.pkl','wb'))
open('progress.log','a').write('PG6 s_C5: delta finalized both Lam conventions; nonzero=%s\n'%(len(d1)>0))
