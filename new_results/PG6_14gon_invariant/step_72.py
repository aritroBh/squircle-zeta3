import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import sympy as sp
g=S.init()
setup=pickle.load(open('q3setup.pkl','rb'))
LI3v=setup['LI3v']; Dgv=setup['Dgv']; basis3=setup['basis3']; CRs=setup['CRs']
bidx={t:i for i,t in enumerate(basis3)}
Mcols=LI3v+Dgv
M=sp.Matrix(27,len(Mcols),lambda i,j: sp.Rational(Mcols[j][i]))
sigma=pickle.load(open('sigma14.pkl','rb'))
s_comp=[[Fr(0)]*27 for _ in range(3)]
for key,cFT in sigma.items():
    A,B,C,d=cFT
    if key in bidx:
        s_comp[0][bidx[key]]+=Fr(A)/d; s_comp[1][bidx[key]]+=Fr(B)/d; s_comp[2][bidx[key]]+=Fr(C)/d
# particular solutions: set free params to 0
sol_c=[[Fr(0)]*24 for _ in range(3)]
for k in range(3):
    tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[k][i]))
    sol,params=M.gauss_jordan_solve(tgt, freevar=True)
    sub={p:sp.Integer(0) for p in params}
    for z in range(24):
        sol_c[k][z]=Fr(sp.Rational(sol[z].subs(sub)) if sol[z].free_symbols else sp.Rational(sol[z]))
print("particular solutions obtained")
# relation nullspace R_li3: kernel of [LI3v | Dgv] restricted to li3 part
# Solve Mli3_rel: columns=LI3v; find c with Mc in Decomp. Build combined and get nullspace of map to Q3.
# Q3 = coker(Decomp). Map phi: Q^24 -> Q3. Kernel = relations.
# Compute: for the 24 li3, phi(c)=0 iff Mc in span(Decomp).
MD = sp.Matrix(27,len(Dgv), lambda i,j: sp.Rational(Dgv[j][i]))
# For each li3 vector, reduce mod Decomp to Q3 coords. Q3 dim 8. Get a basis of Q3 via rref of [LI3|Decomp].
# kernel of phi: solve [LI3v | Dgv] * (c; d) = 0 -> c in relations.
Mnull = M.nullspace()
print("nullspace dim of [li3|Decomp]:", len(Mnull))
# Extract li3-part relations (first 24 components)
Rli3=[]
for vec in Mnull:
    c=[Fr(sp.Rational(vec[z])) for z in range(24)]
    if any(x!=0 for x in c): Rli3.append(c)
print("li3 relations:", len(Rli3))
# delta_li3 for each CR: (vec(z) ^ vec(1-z)) x vec(z)
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
        for l3,c3 in vz.items():
            out[(ab,l3)]=out.get((ab,l3),Fr(0))+cab*c3
    return {k:c for k,c in out.items() if c}
DL=[delta_li3(z) for z in CRs]
def dadd(*ds):
    o={}
    for d in ds:
        for k,c in d.items(): o[k]=o.get(k,Fr(0))+c
    return {k:c for k,c in o.items() if c}
def dsc(d,q):
    q=Fr(q); return {k:c*q for k,c in d.items() if c*q}
# check delta vanishes on relations
def delta_of_coefs(c):
    acc={}
    for z,cz in enumerate(c):
        if cz: acc=dadd(acc, dsc(DL[z],cz))
    return acc
allzero=True
for r in Rli3:
    if delta_of_coefs(r): allzero=False; print("NONZERO delta on relation:", delta_of_coefs(r))
print("delta vanishes on all li3 relations:", allzero)
# delta(sigma) from particular solution
delta_sigma={}
for k in range(3):
    dk=delta_of_coefs(sol_c[k])
    for key,c in dk.items():
        delta_sigma[key]=E.fadd(delta_sigma.get(key,E.ZERO), (Fr(0),)*0 + ( (c,Fr(0),Fr(0),1) if k==0 else (Fr(0),c,Fr(0),1) if k==1 else (Fr(0),Fr(0),c,1)))
print("delta(sigma) terms:", len(delta_sigma))
for key in sorted(delta_sigma): print("  ",key, E.fev(delta_sigma[key],20), delta_sigma[key][:3],"/",delta_sigma[key][3])
pickle.dump({'sol_c':sol_c,'delta_sigma':delta_sigma,'Rli3':Rli3,'DL':DL,'allzero':allzero}, open('delta14.pkl','wb'))
open('progress.log','a').write('PG6 s_B15: delta(sigma) computed, %d terms; delta vanishes on relations: %s\n'%(len(delta_sigma),allzero))
