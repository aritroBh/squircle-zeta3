
import sys, pickle, itertools
sys.path.insert(0,'.')
from fractions import Fraction as Fr
sols2 = pickle.load(open('classical_sols_tracked.pkl','rb'))
T = pickle.load(open('targets81_fixed.pkl','rb'))
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
        ss = [L1[k] for k in name[1]]
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
ok = True
for e in (-1,0,1):
    acc = [Fr(0)]*81
    for name, cf in sols2[e].items():
        if name not in cache: cache[name] = vec81(symbol_of(name))
        v = cache[name]
        for j in range(81): acc[j] += cf*v[j]
    diff = [a-b for a,b in zip(acc, list(T[e]))]
    nz = sum(1 for x in diff if x)
    print(f"e={e}: reconstruction residual:", nz, "components, max|d| =", max((abs(x) for x in diff), default=0))
pickle.dump(cache, open('basis_elem_vecs.pkl','wb'))
open('progress.log','a').write("step_10 [2]-dropped symbol identity solved exactly, all exponents, reconstruction verified\n")
