
import sys, pickle, itertools
sys.path.insert(0,'.')
from fractions import Fraction as Fr
sols = pickle.load(open('classical_sols81_fixed.pkl','rb'))
br = pickle.load(open('basis_red81.pkl','rb'))['basis_red']
T = pickle.load(open('targets81_fixed.pkl','rb'))
# reconstruct: need full vectors for named basis elements -> rebuild from the same construction as before
# basis_red rows are the REDUCED vectors; but solution coefficients refer to the ORIGINAL named basis elements.
# The coefs from solve_target are coefficients of the pivot basis ELEMENTS (by name). Reconstruct by
# re-vectorizing those named elements:
import itertools
seen = pickle.load(open('li_arg_pairs.pkl','rb'))
pairs = sorted(seen.keys())
def shuffle(symA, symB):
    out = {}
    def rec(pref, A, B, cur):
        if not A and not B: out[pref] = out.get(pref, Fr(0)) + cur; return
        if A: rec(pref+(A[0],), A[1:], B, cur)
        if B: rec(pref+(B[0],), A, B[1:], cur)
    for ka, ca in symA.items():
        for kb, cb in symB.items():
            rec((), list(ka), list(kb), ca*cb)
    return {k:v for k,v in out.items() if v}
L1, L2, L3, L4 = {}, {}, {}, {}
for (u, x) in pairs:
    L1[(u,)] = {(u,): Fr(-1)}
    L2[(u,x)] = {(u,x): Fr(-1)}
    L3[(u,x)] = {(u,x,x): Fr(-1)}
    L4[(u,x)] = {(u,x,x,x): Fr(-1)}
def symbol_of(name):
    kind = name[0]
    if kind=='Li4': return L4[name[1]]
    if kind=='L3L1': return shuffle(L3[name[1]], L1[name[2]])
    if kind=='L2L2': return shuffle(L2[name[1]], L2[name[2]])
    if kind=='L2L1L1': return shuffle(shuffle(L2[name[1]], L1[name[2]]), L1[name[3]])
    if kind=='L1^4':
        s = L1[(name[1][0],)] if False else None
        # name[1] is tuple of keys (u,)
        ss = [L1[(k,)] for k in name[1]]
        s = ss[0]
        for j in range(1,4): s = shuffle(s, ss[j])
        return s
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
        v[k[0]*27 + k[1]*9 + k[2]*3 + k[3]] = c0
    return v
cache = {}
for e in (-1,0,1):
    acc = [Fr(0)]*81
    for name, cf in sols[e].items():
        if name not in cache: cache[name] = vec81(symbol_of(name))
        v = cache[name]
        for j in range(81): acc[j] += cf*v[j]
    diff = [a-b for a,b in zip(acc, list(T[e]))]
    print(f"e={e}: reconstruction residual:", sum(1 for x in diff if x), "components, max|d| =", max((abs(x) for x in diff), default=0))
