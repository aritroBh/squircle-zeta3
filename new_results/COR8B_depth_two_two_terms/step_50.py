
import sys, pickle, itertools, time
sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr

seen = pickle.load(open('li_arg_pairs.pkl','rb'))
pairs = sorted(seen.keys())
def shuffle(symA, symB):
    out = {}
    def rec(pref, A, B, cur):
        if not A and not B:
            out[pref] = out.get(pref, Fr(0)) + cur; return
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
basis = {}
for k,s in L4.items(): basis[('Li4',k)] = s
for k1,s1 in L3.items():
    for k2,s2 in L1.items(): basis[('L3L1',k1,k2)] = shuffle(s1,s2)
for i,(k1,s1) in enumerate(L2.items()):
    for k2,s2 in list(L2.items())[i:]: basis[('L2L2',k1,k2)] = shuffle(s1,s2)
for k1,s1 in L2.items():
    for i,(k2,s2) in enumerate(L1.items()):
        for k3,s3 in list(L1.items())[i:]: basis[('L2L1L1',k1,k2,k3)] = shuffle(shuffle(s1,s2),s3)
l1list = list(L1.items())
for combo in itertools.combinations_with_replacement(range(len(l1list)), 4):
    s = l1list[combo[0]][1]
    for j in combo[1:]: s = shuffle(s, l1list[j][1])
    basis[('L1^4', tuple(l1list[j][0] for j in combo))] = s
print("basis size:", len(basis))

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

# incremental RREF (rows as dict-of-dict for sparsity)
t0=time.time()
piv_rows = {}   # pivot_col -> (row dict, name)
def reduce(v):
    v = {i: x for i,x in enumerate(v) if x}
    while v:
        i = min(v)
        if i not in piv_rows:
            return i, v
        row, _ = piv_rows[i]
        f = v[i]/row[i]
        for j,x in row.items():
            if j==i: continue
            nv = v.get(j, Fr(0)) - f*x
            if nv: v[j] = nv
            elif j in v: del v[j]
        del v[i]
    return None, None
basis_red = []  # (name, reduced row dict)
for name, s in basis.items():
    v = vec81(s)
    p, vr = reduce(v)
    if p is not None:
        piv_rows[p] = (vr, name)
        basis_red.append((name, vr))
print("reduced basis (span dim):", len(basis_red), "time", time.time()-t0)
pickle.dump({'basis_red':[(n, dict(r)) for n,r in basis_red]}, open('basis_red81.pkl','wb'))

T = pickle.load(open('targets81.pkl','rb'))
def solve_target(tv):
    # reduce target against piv_rows tracking coefficients
    v = {i: x for i,x in enumerate(tv) if x}
    coef = {}
    while v:
        i = min(v)
        if i not in piv_rows: return None
        row, name = piv_rows[i]
        f = v[i]/row[i]
        coef[name] = coef.get(name, Fr(0)) + f
        for j,x in row.items():
            if j==i: continue
            nv = v.get(j, Fr(0)) - f*x
            if nv: v[j] = nv
            elif j in v: del v[j]
        del v[i]
    return coef
sols = {}
for e in (-1,0,1):
    sol = solve_target(list(T[e]))
    sols[e] = sol
    print("e=",e, "solvable:", sol is not None, "terms:", len(sol) if sol else 0)
pickle.dump(sols, open('classical_sols81.pkl','wb'))
