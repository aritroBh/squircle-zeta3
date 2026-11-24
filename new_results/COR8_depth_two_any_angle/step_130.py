
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
# dedupe to distinct (value, one name)
seen = {}
cands = []
for nm, d in rows:
    key = tuple(sorted((k, tuple(sorted(p.items()))) for k, p in d.items()))
    if key not in seen:
        seen[key] = nm
        cands.append((nm, d))
print(len(cands), "distinct")
target = {((1,2),(2,3)): {0: Fr(1), -1: Fr(-1)},
          ((1,3),(2,3)): {0: Fr(-1), -1: Fr(-1)}}
comps = set()
for nm, d in cands:
    for k, poly in d.items():
        for e in poly:
            comps.add((k, e))
comps = sorted(comps, key=str)
unk = [(j, e) for j in range(len(cands)) for e in (-1, 0, 1)]
A = []
b = []
for (k, e) in comps:
    A.append([sp.Rational(cands[j][1].get(k, {}).get(ee, Fr(0))) for (j, ee) in unk])
    b.append(sp.Rational(target.get(k, {}).get(e, Fr(0))))
M = sp.Matrix(A); bv = sp.Matrix(b)
r = M.row_join(bv).rref()[0]
n_unk = len(unk)
piv_cols = []
for i in range(len(comps)):
    for j in range(n_unk):
        if r[i, j] != 0:
            piv_cols.append(j); break
inconsistent = any(all(r[i, j] == 0 for j in range(n_unk)) and r[i, n_unk] != 0 for i in range(len(comps)))
print("inconsistent:", inconsistent, " pivots:", len(piv_cols))
if not inconsistent:
    sol = [Fr(0)] * n_unk
    for i, j in enumerate(piv_cols):
        sol[j] = Fr(r[i, n_unk])
    used = [(cands[j][0], unk[j][1], sol[j]) for j in range(n_unk) if sol[j]]
    print(len(used), "terms:")
    for nm, e, v in used:
        print(f"  {v} c^{e}  {nm}")
