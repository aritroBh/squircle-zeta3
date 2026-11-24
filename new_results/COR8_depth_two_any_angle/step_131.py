
import pickle
from fractions import Fraction as Fr
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
print("  A=(2m)(cm) B=(2m)(cp) C=(2p)(cm) D=(2p)(cp) E=(cm)(cp) F=(cm)(mp)=W1 G=(cp)(mp)=W2")
for nm, d in cands:
    vec = []
    for k in comps:
        poly = d.get(k, {})
        vec.append('[' + ','.join(str(poly.get(e, '')) for e in (-1,0,1)) + ']')
    print(f"  {nm:44s} " + ' '.join(f"{v:10s}" for v in vec))
