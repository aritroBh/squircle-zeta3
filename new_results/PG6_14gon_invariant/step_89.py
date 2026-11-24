import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import sympy as sp
g=S.init()
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
            if nm.startswith('Lam'):
                sym=lamfun(E.X[int(nm[3])], E.fdiv(E.X[int(nm[5])],E.X[int(nm[3])])); scalar=E.ONE
            else:
                sym,scalar=atom_sym_nonlam(nm)
            tot=E.fmul(c,scalar)
            for key,cf in sym.items():
                sigma[key]=E.fadd(sigma.get(key,E.ZERO), E.fs(tot,cf))
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
rankM=M.rank()
for nm,lamfun in [('mine(-tau)',lam_mine),('dod3-unconv',lam_dod3)]:
    sig=build_sigma(lamfun)
    results=[]
    for k in range(3):
        sc=[Fr(0)]*27
        for key,cFT in sig.items():
            A,B,C,d=cFT
            if key in bidx: sc[bidx[key]]+= [Fr(A)/d,Fr(B)/d,Fr(C)/d][k]
        tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(sc[i]))
        results.append(M.row_join(tgt).rank()==rankM)
    print(f"Lam {nm}: reduces? {results}")
open('progress.log','a').write('PG6 s_C4: Lam convention distinguishing test done\n')
