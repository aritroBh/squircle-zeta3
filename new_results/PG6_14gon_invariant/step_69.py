import pickle, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sym14 as S
import engine14 as E
g = S.init()
# 24 CRs from letters14
L = pickle.load(open('letters14.pkl','rb'))
CRset = L['CRset']
print("CRset entries:", len(CRset), " sample:", CRset[0])
# Each entry is (z, 1-z) as FT. Extract z.
CRs = [e[0] for e in CRset]
# li3 symbols of the 24 CRs (PGEN1R convention z x z x (1-z))
def li3sym(z):
    return S.li3_symbol(z)
LI3 = [li3sym(z) for z in CRs]
# quotient: V={0,1,2}, basis3=27, Decomp
letters=[0,1,2]
basis3=list(itertools.product(letters,repeat=3)); bidx={t:i for i,t in enumerate(basis3)}
Dg=[]
for a in letters:
    for b in letters:
        for c in letters:
            v={}
            for t in [(a,b,c),(b,a,c),(b,c,a)]:
                v[t]=v.get(t,0)+1
            Dg.append(v)
# represent tensors as 27-dim rational vectors
def t2v(T):
    v=[Fr(0)]*27
    for t,c in T.items():
        if t in bidx: v[bidx[t]]+=Fr(c)
    return v
import sympy as sp
DM = sp.Matrix(27, len(Dg), lambda i,j: sp.Rational(Dg[j].get(basis3[i],0)))
print("Decomp rank:", DM.rank(), " (expect 19), dim Q3 =", 27-DM.rank())
# build the 24 li3 vectors
LI3v=[t2v(T) for T in LI3]
# reduce sigma (FT coeffs) -> 3 rational components
sigma = pickle.load(open('sigma14.pkl','rb'))
# sigma: key=(i,j,k), coef FT. Convert to 3 rational 27-vectors (1,alpha,alpha^2 comps)
s_comp=[[Fr(0)]*27 for _ in range(3)]
for key,cFT in sigma.items():
    A,B,C,d = cFT
    if key in bidx:
        s_comp[0][bidx[key]] += Fr(A)/d
        s_comp[1][bidx[key]] += Fr(B)/d
        s_comp[2][bidx[key]] += Fr(C)/d
# Check sigma is in span(li3)+Decomp, solve over Q for each component
Mcols = LI3v + Dg  # 24 + 27 columns
M = sp.Matrix(27, len(Mcols), lambda i,j: sp.Rational(Mcols[j][i]))
print("rank [li3|Decomp]:", M.rank())
pickle.dump({'LI3v':LI3v,'Dg':Dg,'basis3':basis3,'CRs':CRs}, open('q3setup.pkl','wb'))
open('progress.log','a').write('PG6 s_B12: Q3 quotient built (dim 8), 24 CR li3 symbols; rank[li3|Decomp]=%d\n'%M.rank())
