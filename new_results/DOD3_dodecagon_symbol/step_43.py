import sys, os, pickle, itertools
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
import sympy as sp
from dod3_symbol import *

def shuffle3(u, v, w):
    vecs = [u,v,w]
    S = {}
    for perm in itertools.permutations(range(3)):
        U,V,W = vecs[perm[0]], vecs[perm[1]], vecs[perm[2]]
        for a,ca in U.items():
            for b,cb in V.items():
                for c,cc in W.items():
                    S[(a,b,c)] = S.get((a,b,c),Fr(0)) + ca*cb*cc
    return {k:x for k,x in S.items() if x!=0}

def S_logprod(letters):
    vecs = [{l: Fr(1)} for l in letters]
    if len(vecs)==1:
        u = vecs[0]; return tensor3(u,u,u)
    if len(vecs)==2:
        u,v = vecs
        S = {}
        for a,ca in u.items():
            for b,cb in v.items():
                S[(a,b,b)] = S.get((a,b,b),Fr(0))+ca*cb
                S[(b,a,b)] = S.get((b,a,b),Fr(0))+ca*cb
                S[(b,b,a)] = S.get((b,b,a),Fr(0))+ca*cb
        return {k:x for k,x in S.items() if x!=0}
    return shuffle3(*vecs)

def S_Li2log(w, p):
    wv = vecfun(w); owv = vecfun(psub(ONE,w))
    S2 = {}
    for a,ca in owv.items():
        for b,cb in wv.items():
            S2[(a,b)] = S2.get((a,b),Fr(0)) - ca*cb
    S = {}
    for (l1,l2_),c in S2.items():
        for t in [(p,l1,l2_),(l1,p,l2_),(l1,l2_,p)]:
            S[t] = S.get(t,Fr(0)) + c
    return {k:x for k,x in S.items() if x!=0}

# sanity: S_logprod(['2','3']) should equal symbol of l2*l3 weight-2: 2t3 + 3t2
print("S_logprod(['2','3']):", S_logprod(['2','3']))
prod_syms = []; prod_names = []
for a in range(4):
    for b in range(4-a):
        c = 3-a-b
        letters = ['2']*a + ['3']*b + ['e']*c
        prod_syms.append(S_logprod(letters))
        prod_names.append("l2^%d l3^%d lu^%d"%(a,b,c))
b2zs = pickle.load(open('li2_basis.pkl','rb'))['b2_zs']
for w in b2zs:
    for p in ('2','3','e'):
        prod_syms.append(S_Li2log(w,p))
        prod_names.append("Li2(%s)*l%s"%(str(w),p))
Q = Quotient()
Pcols = [Q.coord(s) for s in prod_syms]
PM = sp.Matrix.hstack(*[c.T for c in Pcols])
print("rank of 19 product symbols:", PM.rank(), "(expect 19)")
# check span equals Decomp: rank [PM | D] should be 19
DM_full = sp.Matrix.hstack(*[d.T for d in Q.D])
print("rank [products | Decomp]:", sp.Matrix.hstack(PM, DM_full).rank(), "(expect 19)")
pickle.dump({'prod_syms':prod_syms,'prod_names':prod_names}, open('prod_sym_data.pkl','wb'))
