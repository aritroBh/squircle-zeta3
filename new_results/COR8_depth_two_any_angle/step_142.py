
import sys, itertools, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import itertools
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/symbols_2kept.pkl','rb'))
OM, DT, R = D['OM'], D['DT'], D['R']
bwords = list(itertools.product(range(4), repeat=4))
# does OM have [2]-content?
om2 = [(w, cf) for w, cf in OM.items() if any(l[0] for l in w)]
print("OM words with [2]:", len(om2))
# DT words with [2]
def ws(w):
    out = []
    for l in w:
        if sum(map(abs, l)) == 1:
            out.append('2cmp'[next(j for j, x in enumerate(l) if x)])
        else:
            out.append(str(tuple(int(x) for x in l)))
    return ''.join(out)
dt2 = [(w, cf) for w, cf in DT.items() if any(l[0] for l in w)]
print("DT words with [2]:", len(dt2))
for w, cf in sorted(dt2, key=str)[:10]:
    print("  ", ws(w), {e: str(v) for e, v in sorted(cf.items())})
# check R = OM - DT at a [2] word
for probe in dt2[:3]:
    w, cf = probe
    print("R at", ws(w), {e: str(v) for e, v in sorted(R.get(w, {}).items())})
