
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
# keep only rows whose delta22 has some [2] content
rows2 = [(nm, d) for nm, d in rows if any(0 in (kk[0]+kk[1]) for kk in d)]
print(len(rows2), "pairs with [2]-content in delta22")
target = {((0,2),(1,2)): {0: Fr(1), -1: Fr(-1)},
          ((0,2),(1,3)): {0: Fr(1), -1: Fr(1)},
          ((0,3),(1,2)): {0: Fr(1), -1: Fr(-1)},
          ((0,3),(1,3)): {0: Fr(1), -1: Fr(1)}}
# solve sum r_j d_j = target with r_j Laurent (e in {-1,0,1}), non-target comps = 0
comps = set()
for nm, d in rows2:
    for k, poly in d.items():
        for e in poly:
            comps.add((k, e))
comps = sorted(comps, key=str)
unk = [(j, e) for j in range(len(rows2)) for e in (-1, 0, 1)]
# equations: for comp in comps: sum = target
import sympy as sp
A = []
b = []
for (k, e) in comps:
    row = []
    for (j, ee) in unk:
        row.append(sp.Rational(rows2[j][1].get(k, {}).get(ee, Fr(0))))
    A.append(row)
    b.append(sp.Rational(target.get(k, {}).get(e, Fr(0))))
M = sp.Matrix(A)
bv = sp.Matrix(b)
aug = M.row_join(bv)
r = aug.rref()[0]
# check consistency
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
    used = [(rows2[j][0], unk[j][1], sol[j]) for j in range(n_unk) if sol[j]]
    print(len(used), "terms:")
    for nm, e, v in used:
        print(f"  {v} c^{e}  {nm}")
