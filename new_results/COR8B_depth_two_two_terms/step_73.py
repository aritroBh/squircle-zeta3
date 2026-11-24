
import sys, pickle, itertools, time
sys.path.insert(0,'.')
from fractions import Fraction as Fr
S = pickle.load(open('syms_omega_dt.pkl','rb'))
S_diff = S['S_diff']
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
def target_vec_fixed(e):
    out = {}
    for k, poly in S_diff.items():
        c0 = poly.get(e, Fr(0))
        if not c0: continue
        nk = []
        ok = True
        for letter in k:
            p = letter[1:]          # PROJECT: set [2]-component to zero
            if all(x==0 for x in p): ok = False; break
            nk.append(p)
        if not ok: continue
        nk = tuple(nk)
        out[nk] = out.get(nk, Fr(0)) + c0
    return vec81(out)
T = {e: target_vec_fixed(e) for e in (-1,0,1)}
print("fixed target component counts:", {e: sum(1 for x in v if x) for e,v in T.items()})
pickle.dump(T, open('targets81_fixed.pkl','wb'))

br = pickle.load(open('basis_red81.pkl','rb'))['basis_red']
piv_rows = {}
for name, row in br:
    piv_rows[min(row)] = (row, name)
def solve_target(tv):
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
    print("e=",e,"solvable:", sol is not None, "terms:", len(sol) if sol else 0)
pickle.dump(sols, open('classical_sols81_fixed.pkl','wb'))
