import pickle, itertools
from fractions import Fraction as Fr
import sympy as sp
from q1s_symbol import tadd, t_conj
from q1s_exact import fconj, feq
B4=('a','p2','q1','q2')
BIDX=list(itertools.product(B4,repeat=3))
def coord(T):
    v=sp.zeros(64,1)
    for k,c in T.items(): v[BIDX.index(k)]=sp.Rational(c)
    return v
def deco_gens():
    gens=[]
    for u in B4:
        for v in B4:
            for w in B4:
                t={}
                for perm in [(u,v,w),(v,u,w),(v,w,u)]: t[perm]=t.get(perm,Fr(0))+1
                gens.append(coord(t))
    return gens
D=deco_gens()
sols=pickle.load(open('sunit_li3.pkl','rb'))
pairs=[]; used=set()
for i,s in enumerate(sols):
    if i in used: continue
    j=next((jj for jj,s2 in enumerate(sols) if jj!=i and feq(s2['z'],fconj(s['z']))), None)
    used.add(i)
    if j is not None: used.add(j)
    pairs.append((i,j))
Gp=[]
for i,j in pairs:
    g=tadd(sols[i]['S'], t_conj(sols[i]['S']))
    Gp.append((i,j,coord(g)))
Gr=[g for _,_,g in Gp]
Mr=sp.Matrix.hstack(*Gr)
print("rank[272 ReLi3]:", Mr.rank())
Mrd=sp.Matrix.hstack(*(Gr+D))
print("rank[ReLi3|Decomp]:", Mrd.rank())
S=pickle.load(open('q3_symbols.pkl','rb'))
for nm,t in [("V3'",coord(S['S3_V3'][0])),("V1 A",coord(S['S3_V1'][0])),("V1 B",coord(S['S3_V1'][1]))]:
    print(nm, "in real span:", sp.Matrix.hstack(Mrd,t).rank()==Mrd.rank())
