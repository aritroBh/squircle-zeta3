import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
import sympy as sp
g = S.init()
L = pickle.load(open('letters14.pkl','rb'))
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
print("rank[li3|Decomp]:", M.rank())
sigma=pickle.load(open('sigma14.pkl','rb'))
s_comp=[[Fr(0)]*27 for _ in range(3)]
for key,cFT in sigma.items():
    A,B,C,d=cFT
    if key in bidx:
        s_comp[0][bidx[key]]+=Fr(A)/d; s_comp[1][bidx[key]]+=Fr(B)/d; s_comp[2][bidx[key]]+=Fr(C)/d
rankM=M.rank(); sol_c=[[None]*24 for _ in range(3)]; ok=True
for k in range(3):
    tgt=sp.Matrix(27,1,lambda i,j: sp.Rational(s_comp[k][i]))
    if M.row_join(tgt).rank()!=rankM:
        print(f"comp {k}: NOT in span"); ok=False; continue
    sol=M.gauss_jordan_solve(tgt)[0]
    for z in range(24): sol_c[k][z]=Fr(sp.Rational(sol[z]))
print("all comps in span:", ok)
pickle.dump({'LI3v':LI3v,'Dgv':Dgv,'basis3':basis3,'CRs':CRs,'sol_c':sol_c,'ok':ok}, open('q3setup.pkl','wb'))
# print the c_z
if ok:
    print("c_z (as FT) for the 24 CRs:")
    for z in range(24):
        c=(sol_c[0][z],sol_c[1][z],sol_c[2][z])
        if any(x!=0 for x in c):
            print("  z=%d"%z, c, " eval", E.fev((c[0],c[1],c[2],1),20))
open('progress.log','a').write('PG6 s_B14: sigma reduced to li3 basis over F: ok=%s\n'%ok)
