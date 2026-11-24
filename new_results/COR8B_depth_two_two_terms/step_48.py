
import sys, pickle, itertools, time
sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr

S = pickle.load(open('syms_omega_dt.pkl','rb'))
S_diff = S['S_diff']
seen = pickle.load(open('li_arg_pairs.pkl','rb'))
pairs = sorted(seen.keys())   # ([1-g],[g]) as 3-vectors (c,p,q coords)

def T3(*vs): return {tuple(vs): Fr(1)}
def shuffle(symA, symB):
    # recursive shuffle of two symbols (any weights)
    la = len(next(iter(symA))) if symA else 0
    lb = len(next(iter(symB))) if symB else 0
    out = {}
    def rec(pref, A, B):
        if not A and not B:
            out[pref] = out.get(pref, Fr(0)) + cur
            return
        if A:
            rec(pref+(A[0],), A[1:], B)
        if B:
            rec(pref+(B[0],), A, B[1:])
    for ka, ca in symA.items():
        for kb, cb in symB.items():
            cur = ca*cb
            rec((), list(ka), list(kb))
    return {k:v for k,v in out.items() if v}

L1, L2, L3, L4 = {}, {}, {}, {}
for (u, x) in pairs:  # u=[1-g], x=[g]
    L1[(u,)] = {(u,): Fr(-1)}           # S(Li1(g)) = -[1-g]
    L2[(u,x)] = {(u,x): Fr(-1)}         # S(Li2(g)) = -[1-g][g]
    L3[(u,x)] = {(u,x,x): Fr(-1)}       # S(Li3(g))
    L4[(u,x)] = {(u,x,x,x): Fr(-1)}     # S(Li4(g))
print("L1,L2,L3,L4:", len(L1),len(L2),len(L3),len(L4))

# product basis
basis = {}   # name -> symbol
for k,s in L4.items(): basis[('Li4',k)] = s
t0=time.time()
for k1,s1 in L3.items():
    for k2,s2 in L1.items():
        basis[('Li3*Li1',k1,k2)] = shuffle(s1,s2)
for i,(k1,s1) in enumerate(L2.items()):
    for k2,s2 in list(L2.items())[i:]:
        basis[('Li2*Li2',k1,k2)] = shuffle(s1,s2)
for k1,s1 in L2.items():
    for i,(k2,s2) in enumerate(L1.items()):
        for k3,s3 in list(L1.items())[i:]:
            basis[('Li2*Li1*Li1',k1,k2,k3)] = shuffle(shuffle(s1,s2),s3)
l1list = list(L1.items())
for combo in itertools.combinations_with_replacement(range(len(l1list)), 4):
    s = l1list[combo[0]][1]
    for j in combo[1:]:
        s = shuffle(s, l1list[j][1])
    basis[('Li1^4', tuple(l1list[j][0] for j in combo))] = s
print("basis size:", len(basis), "time", time.time()-t0)

# vector space: all length-4 tensors over 3 letters -> 81-dim
def to_vec(sym):
    v = {}
    for k,c0 in sym.items():
        idx = 0
        for letter in k:
            # letter is 3-vector; but our symbols are over arbitrary 3-vectors -> treat each distinct letter as coordinate? 
            pass
    return v
# letters are arbitrary integer 3-vectors; the tensor algebra is over the 3-dim space span{c,p,q}
# -> convert: each letter v=(a,b,d) = a*e1+b*e2+d*e3, expand multilinearly
E = [(1,0,0),(0,1,0),(0,0,1)]
def expand(sym):
    out = {}
    for k,c0 in sym.items():
        opts = [[(i,vi) for i,vi in enumerate(letter) if vi] for letter in k]
        for combo in itertools.product(*opts):
            mult = Fr(c0)
            for i,vi in combo: mult *= vi
            nk = tuple(i for i,vi in combo)
            out[nk] = out.get(nk, Fr(0)) + mult
        out = {k:v for k,v in out.items() if v}
    return out
def vec81(sym):
    e = expand(sym)
    v = [Fr(0)]*81
    for k,c0 in e.items():
        idx = k[0]*27 + k[1]*9 + k[2]*3 + k[3]
        v[idx] = c0
    return tuple(v)
Bvec = {name: vec81(s) for name,s in basis.items()}
print("basis vectorized")
# target vectors per exponent
def target_vec(e):
    out = {}
    for k, poly in S_diff.items():
        if any(letter[0] != 0 for letter in k): continue  # drop [2]-tensors
        c0 = poly.get(e, Fr(0))
        if c0: out[k] = out.get(k, Fr(0)) + c0
    return vec81(out)
T = {e: target_vec(e) for e in (-1,0,1)}
print("target component counts:", {e: sum(1 for x in v if x) for e,v in T.items()})
pickle.dump({'Bvec':Bvec,'T':T}, open('solve81_data.pkl','wb'))
# rank check
import sympy
M = sympy.Matrix([list(v) for v in Bvec.values()])
t0=time.time(); r = M.rank(); print("basis rank:", r, "time", time.time()-t0)
for e in (-1,0,1):
    Mt = sympy.Matrix([list(v) for v in Bvec.values()] + [list(T[e])])
    print("e=",e,"target in span:", Mt.rank()==r)
