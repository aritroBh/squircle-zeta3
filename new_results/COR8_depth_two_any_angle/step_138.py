
import pickle
from fractions import Fraction as Fr
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/symbols_2kept.pkl','rb'))
R = D['R']
def ws(w):
    return ''.join('2cmp'[next(j for j,x in enumerate(l) if x)] if sum(map(abs,l))==1 else str(tuple(int(x) for x in l)) for l in w)
with2 = [(w, cf) for w, cf in R.items() if any(l[0] for l in w)]
print(len(with2), "words of R containing [2]")
for w, cf in sorted(with2, key=lambda x: str(x[0])):
    print("  ", ws(w), {e: str(v) for e, v in sorted(cf.items())})
