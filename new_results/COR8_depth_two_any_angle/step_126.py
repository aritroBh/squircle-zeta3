
import pickle
from fractions import Fraction as Fr
import sympy as sp
rows = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/pairs_delta22_2kept.pkl','rb'))
def parse_d(d):
    out = {}
    for k, poly in d.items():
        kk = eval(k)
        out[kk] = {eval(e): Fr(v) for e, v in poly.items()}
    return out
rows = [(nm, parse_d(d)) for nm, d in rows]
target_keys = {((0,2),(1,2)), ((0,2),(1,3)), ((0,3),(1,2)), ((0,3),(1,3))}
sub = [(nm, d) for nm, d in rows if set(d.keys()) <= target_keys]
print(len(sub), "pairs with support inside target")
for nm, d in sub:
    print(nm, {k: {e: str(v) for e,v in sorted(p.items())} for k,p in sorted(d.items())})
