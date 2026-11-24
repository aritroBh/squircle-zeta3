
import sys, pickle, itertools, time
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import sympy as sp
c = sp.symbols('c')
vals = pickle.load(open('sunit_args_fixed.pkl','rb'))
SAFE = ['0.500000000000000','0.5*c + 0.5','0.5 - 0.5*c','1/(c + 1)','(1 - c)/(c + 1)','-(1 - c)/(c + 1)',
        '1 - c','1 - c**2','(c**2 - 2*c + 1)/(c**2 + 2*c + 1)','c/(c + 1)','c','-c','c**2','2*c/(c + 1)',
        '4*c/(c**2 + 2*c + 1)','1','-1']
safe_pairs = {}
for a in SAFE:
    vt, v1t = vals[a]
    safe_pairs.setdefault((v1t[1:], vt[1:]), sp.sympify(a))
print("safe distinct pairs:", len(safe_pairs))
pairs = sorted(safe_pairs.keys())
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
print("L1,L2,L3,L4:", len(L1),len(L2),len(L3),len(L4))
basis = []
for k,s in L4.items(): basis.append((('Li4',k), s))
for k1,s1 in L3.items():
    for k2,s2 in L1.items(): basis.append((('L3L1',k1,k2), shuffle(s1,s2)))
for i,(k1,s1) in enumerate(L2.items()):
    for k2,s2 in list(L2.items())[i:]: basis.append((('L2L2',k1,k2), shuffle(s1,s2)))
for k1,s1 in L2.items():
    for i,(k2,s2) in enumerate(L1.items()):
        for k3,s3 in list(L1.items())[i:]: basis.append((('L2L1L1',k1,k2,k3), shuffle(shuffle(s1,s2),s3)))
l1list = list(L1.items())
for combo in itertools.combinations_with_replacement(range(len(l1list)), 4):
    s = l1list[combo[0]][1]
    for j in combo[1:]: s = shuffle(s, l1list[j][1])
    basis.append((('L1^4', tuple(l1list[j][0] for j in combo)), s))
print("basis:", len(basis))
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
t0=time.time()
piv_rows = {}
def reduce_tracked(v, cdict):
    v = {i: x for i,x in enumerate(v) if x}
    cd = dict(cdict)
    while v:
        i = min(v)
        if i not in piv_rows: return i, v, cd
        row, rc = piv_rows[i]
        f = v[i]/row[i]
        for j,x in row.items():
            if j==i: continue
            nv = v.get(j, Fr(0)) - f*x
            if nv: v[j] = nv
            elif j in v: del v[j]
        del v[i]
        for nm,x in rc.items():
            nc = cd.get(nm, Fr(0)) - f*x
            if nc: cd[nm] = nc
            elif nm in cd: del cd[nm]
    return None, None, None
for name, s in basis:
    p, vr, cd = reduce_tracked(vec81(s), {name: Fr(1)})
    if p is not None: piv_rows[p] = (vr, cd)
print("safe span dim:", len(piv_rows), "time", time.time()-t0)
pickle.dump({'piv':{p:(dict(r),dict(cd)) for p,(r,cd) in piv_rows.items()},'basis_meta':'safe'}, open('piv_tracked_safe.pkl','wb'))
T = pickle.load(open('targets81_fixed.pkl','rb'))
def solve(tv):
    v = {i: x for i,x in enumerate(tv) if x}
    cd = {}
    while v:
        i = min(v)
        if i not in piv_rows: return None
        row, rc = piv_rows[i]
        f = v[i]/row[i]
        for j,x in row.items():
            if j==i: continue
            nv = v.get(j, Fr(0)) - f*x
            if nv: v[j] = nv
            elif j in v: del v[j]
        del v[i]
        for nm,x in rc.items():
            cd[nm] = cd.get(nm, Fr(0)) + f*x
    return cd
sols_safe = {}
for e in (-1,0,1):
    sol = solve(list(T[e]))
    sols_safe[e] = sol
    print("e=",e,"solvable:", sol is not None, "terms:", len(sol) if sol else 0)
pickle.dump(sols_safe, open('classical_sols_safe.pkl','wb'))
pickle.dump({str(k): str(v) for k,v in safe_pairs.items()}, open('safe_pairs.pkl','wb'))
