import pickle, itertools
from fractions import Fraction as Fr
import sympy as sp
from q1s_symbol import tadd, tscale, t_conj
import time
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
t0=time.time()
Gc=[coord(s['S']) for s in sols]
Mc=sp.Matrix.hstack(*Gc)
print("rank[493 Li3]:", Mc.rank())
MD=sp.Matrix.hstack(*D)
print("rank[Decomp]:", MD.rank())
Mall=sp.Matrix.hstack(*(Gc+D))
print("rank[all]:", Mall.rank(), " %.1fs"%(time.time()-t0))
# targets
S=pickle.load(open('q3_symbols.pkl','rb'))
tV3=coord(S['S3_V3'][0])
tV1A=coord(S['S3_V1'][0]); tV1B=coord(S['S3_V1'][1])
for nm,t in [("V3'",tV3),("V1 A",tV1A),("V1 B",tV1B)]:
    Mx=sp.Matrix.hstack(Mall, t)
    print(nm, "in complex span:", Mx.rank()==Mall.rank())
