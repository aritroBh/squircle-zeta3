
import sys, itertools, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import itertools
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve256.pkl','rb'))
Rv = D['Rv']
bwords = list(itertools.product(range(4), repeat=4))
bidx = {w: i for i, w in enumerate(bwords)}
def widx(s):
    return bidx[tuple('2cmp'.index(ch) for ch in s)]
func = {'222m': 1, '22m2': -1, '22mm': -1, '22mp': -1, 'm22m': 1, 'p22m': 1}
for E in (-1, 0, 1):
    val = sum(c * Rv[E][widx(w)] for w, c in func.items())
    print(f"f . R_{E} = {val}")
# check content of R_-1 at those words
for w, cf in func.items():
    print(w, Rv[-1][widx(w)])
