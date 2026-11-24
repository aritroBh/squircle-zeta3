
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
M = sp.Matrix([[sp.Rational(cands[j][1].get(k, {}).get(0, Fr(0))) for j in range(len(cands))] for k in comps])
print("rank full 7x32:", M.rank())
# kernel constraints on (F,G) image: left null space of rows A..E against (F,G)
# achievable (F,G) = image of M restricted to ker(A..E). Compute: null space of M[0:5,:] in R^32, then project to (F,G).
ker = M[:5, :].nullspace()
print("dim ker(A..E):", len(ker))
FG = []
for kv in ker:
    f = sum(kv[j]*M[5, j] for j in range(len(cands)))
    g = sum(kv[j]*M[6, j] for j in range(len(cands)))
    FG.append((sp.nsimplify(f), sp.nsimplify(g)))
print("achievable (F,G) with A..E=0:", set(FG))
# row space image basis
img = sp.Matrix.hstack(*[kv for kv in ker]).T
proj = sp.Matrix([[sum(kv[j]*M[5, j] for j in range(len(cands))), sum(kv[j]*M[6, j] for j in range(len(cands)))] for kv in ker])
print("(F,G) image rank:", proj.rank(), "rows:", proj.tolist())
