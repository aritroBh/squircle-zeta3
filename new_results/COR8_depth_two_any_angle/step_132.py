
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
# add the two DT pair values explicitly if missing
from cor8_sym import *
Z4 = (Fr(0),)*4
def v4(j, val=1):
    v = list(Z4); v[j] = Fr(val); return tuple(v)
L2, cv, mv, pv = v4(0), v4(1), v4(2), v4(3)
def dropz(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
Rm = (Fr(-1), Fr(0), Fr(0), Fr(1)); Rp = (Fr(-1), Fr(0), Fr(1), Fr(0))
d22m = delta22(dropz(li22_sym(L2, mv, L2, cv, Rm)), 4)
d22p = delta22(dropz(li22_sym(L2, pv, L2, cv, Rp)), 4)
extra = [('Li22(1/2,1-c)', d22m), ('Li22(1/2,1+c)', d22p)]
for nm, d in extra:
    key = tuple(sorted((k, tuple(sorted(p.items()))) for k, p in d.items()))
    if key not in seen:
        seen[key] = nm; cands.append((nm, d))
print(len(cands), "candidates")
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
tgt_const = {((1,2),(2,3)): Fr(1), ((1,3),(2,3)): Fr(-1)}
tgt_cinv  = {((1,2),(2,3)): Fr(-1), ((1,3),(2,3)): Fr(-1)}
s_const = solve(tgt_const)
s_cinv = solve(tgt_cinv)
print("const solution:", s_const is not None, " 1/c solution:", s_cinv is not None)
if s_const:
    for nm, v in [(cands[j][0], s_const[j]) for j in range(len(cands)) if s_const[j]]:
        print("  ", v, nm)
if s_cinv:
    for nm, v in [(cands[j][0], s_cinv[j]) for j in range(len(cands)) if s_cinv[j]]:
        print("  ", v, "c^-1 ", nm)
