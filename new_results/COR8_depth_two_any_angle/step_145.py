
import pickle
from fractions import Fraction as Fr
import sympy as sp
rows = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/pairs_delta22_2kept.pkl','rb'))
def parse_d(d):
    out = {}
    for k, poly in d.items():
        out[eval(k)] = {eval(e): Fr(v) for e, v in poly.items()}
    return out
rows = [(nm, parse_d(d)) for nm, d in rows]
seen = {}
cands = []
for nm, d in rows:
    key = tuple(sorted((k, tuple(sorted(p.items()))) for k, p in d.items()))
    if key not in seen:
        seen[key] = nm
        cands.append((nm, d))
comps = [((0,2),(1,2)), ((0,2),(1,3)), ((0,3),(1,2)), ((0,3),(1,3)), ((1,2),(1,3)), ((1,2),(2,3)), ((1,3),(2,3))]
def solve(tgt):
    A = sp.Matrix([[sp.Rational(cands[j][1].get(k, {}).get(0, Fr(0))) for j in range(len(cands))] for k in comps])
    b = sp.Matrix([sp.Rational(tgt.get(k, Fr(0))) for k in comps])
    r = A.row_join(b).rref()[0]
    n = len(cands)
    piv = []
    for i in range(len(comps)):
        for j in range(n):
            if r[i,j]: piv.append(j); break
    inc = any(all(r[i,j]==0 for j in range(n)) and r[i,n]!=0 for i in range(len(comps)))
    if inc: return None
    sol = [Fr(0)]*n
    for i, j in enumerate(piv): sol[j] = Fr(r[i, n])
    return sol
# target for extra depth-two: A: -(1+1/c), B: -(1+1/c), C: -(1-1/c), D: (1+1/c)
tgt = {comps[0]: -1, comps[1]: -1, comps[2]: -1, comps[3]: 1}
tgt2 = {comps[0]: -1, comps[1]: -1, comps[2]: 1, comps[3]: 1}
for label, t in [('const', tgt), ('1/c', tgt2)]:
    s = solve(t)
    print(label, "solution:", s is not None)
    if s:
        for nm, v in [(cands[j][0], s[j]) for j in range(len(cands)) if s[j]]:
            print("   ", v, nm)
# also try with F,G allowed to be nonzero (must cancel among extras)
