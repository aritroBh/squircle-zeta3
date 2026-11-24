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
# build sigma as 3 rational 27-vectors
s_comp=[[Fr(0)]*27 for _ in range(3)]
for key,cFT in sigma.items():
    A,B,C,d=cFT
    if key in bidx: s_comp[0][bidx[key]]+=Fr(A)/d; s_comp[1][bidx[key]]+=Fr(B)/d; s_comp[2][bidx[key]]+=Fr(C)/d
# robust solve: for each component, find c_z (24) and decomp (27) with [LI3|D]*(c;d)=target.
# Use particular solution via pseudoinverse-free rref: solve ML3*c + MD*d = tgt. 
# Then VERIFY: rank([ML3c-tgt | MD]) == rank(MD) (i.e., residual in Decomp).
MD=sp.Matrix(27,len(Dgv),lambda i,j: sp.Rational(Dgv[j][i]))
def solve_comp(k, seed=0):
    tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[k][i]))
    # gauss_jordan_solve, extract free params properly, then random-set them
    res=M.gauss_jordan_solve(tgt)
    sol=res[0]
    # free params
    import itertools as it
    free=set()
    for z in range(sol.rows): free|=sol[z].free_symbols
    free=sorted(free, key=str)
    sub={fs:sp.Integer((z*seed+seed)%7-3) for z,fs in enumerate(free)}
    solp=[Fr(sp.Rational(sol[z].subs(sub))) for z in range(24)]
    return solp
def delta_of(c):
    acc={}
    for z,cz in enumerate(c):
        if cz:
            for key,cf in DL[z].items():
                acc[key]=acc.get(key,Fr(0))+cf*cz
    return {k:v for k,v in acc.items() if v}
# compute delta for component 0 with 4 different seeds (particular solutions)
for seed in [0,1,2,3]:
    c0=solve_comp(0,seed)
    d0=delta_of(c0)
    # verify solution: residual in Decomp
    resid = ML3*sp.Matrix(24,1,lambda i,j: sp.Rational(c0[i])) - sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[0][i]))
    rk = MD.row_join(resid).rank()
    print(f"seed {seed}: delta comp0 = {d0}, resid-in-Decomp: {rk==MD.rank()}")
