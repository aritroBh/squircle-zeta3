
import pickle
from fractions import Fraction as Fr
rows = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/pairs_delta22_2kept.pkl','rb'))
def parse_d(d):
    out = {}
    for k, poly in d.items():
        out[eval(k)] = {eval(e): Fr(v) for e, v in poly.items()}
    return out
rows = [(nm, parse_d(d)) for nm, d in rows]
from collections import Counter, defaultdict
supp = Counter()
byval = defaultdict(list)
for nm, d in rows:
    supp[tuple(sorted(d.keys()))] += 1
    key = tuple(sorted((k, tuple(sorted(p.items()))) for k, p in d.items()))
    byval[key].append(nm)
for s, n in supp.most_common():
    print(s, n)
print()
print("distinct delta22 values:", len(byval))
for key, occ in sorted(byval.items(), key=lambda kv: -len(kv[1]))[:14]:
    d = {k: dict(v) for k, v in key}
    print({k: {e: str(v) for e,v in sorted(vv.items())} for k,vv in d.items()}, "x%d" % len(occ), occ[:3])
