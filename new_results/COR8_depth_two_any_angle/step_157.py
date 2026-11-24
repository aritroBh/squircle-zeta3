
import pickle, itertools
from fractions import Fraction as Fr
Rv = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/Rv_2kept.pkl','rb'))
bwords = list(itertools.product(range(4), repeat=4))
bidx = {w: i for i, w in enumerate(bwords)}
def widx(s):
    return bidx[tuple('2cmp'.index(ch) for ch in s)]
F1 = {'m2cm': 1, 'mcm2': -1, 'mcmp': -1, 'mcpm': 1}
F2 = {'2mcm': 1, '2pcm': 1, 'm2cm': -1, 'pcm2': -1, 'pcmp': -1, 'pcpm': 1}
for E in (-1, 0, 1):
    for nm, F in [('F1', F1), ('F2', F2)]:
        val = sum(co*Rv[E][widx(w)] for w, co in F.items())
        if val:
            print(f"E={E} {nm}: {val}")
            for w, co in F.items():
                if Rv[E][widx(w)]:
                    print(f"    {w}: {Rv[E][widx(w)]}")
