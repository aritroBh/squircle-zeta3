import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import sympy as sp
g=S.init()
# quotient setup
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
# full relation space R: c in Q^24 with (sum c li3) in Decomp
# = projection of null([LI3v|Dgv]) onto li3 coords
M=sp.Matrix(27,len(LI3v)+len(Dgv),lambda i,j: sp.Rational((LI3v+Dgv)[j][i]))
ns=M.nullspace()
# relation: for each null vector, li3 part (first 24)
Rli3=[]
for vec in ns:
    c=[Fr(sp.Rational(vec[z])) for z in range(24)]
    if any(x!=0 for x in c): Rli3.append(c)
# reduce Rli3 to independent set (row rank)
if Rli3:
    RM=sp.Matrix(len(Rli3),24,lambda i,j: sp.Rational(Rli3[i][j]))
    print("dim relation space (li3 combos in Decomp):", RM.rank(), " (expect 24 - rank_li3modD)")
ML3=sp.Matrix(27,24,lambda i,j: sp.Rational(LI3v[j][i]))
MD=sp.Matrix(27,len(Dgv),lambda i,j: sp.Rational(Dgv[j][i]))
# rank of li3 mod Decomp = rank([LI3|D]) - rank(D)
print("rank li3 mod Decomp:", M.rank()-MD.rank(), " dim Q3:", 27-MD.rank())
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
def dadd(*ds):
    o={}
    for d in ds:
        for k,c in d.items(): o[k]=o.get(k,Fr(0))+c
    return {k:c for k,c in o.items() if c}
def dsc(d,q):
    q=Fr(q); return {k:c*q for k,c in d.items() if c*q}
def delta_of(c):
    acc={}
    for z,cz in enumerate(c):
        if cz: acc=dadd(acc,dsc(DL[z],cz))
    return acc
# Check delta on EACH independent relation
bad=0
for r in Rli3:
    dr=delta_of(r)
    if dr: bad+=1; print("RELATION WITH NONZERO DELTA:", dr)
print("relations with nonzero delta:", bad, "of", len(Rli3))
