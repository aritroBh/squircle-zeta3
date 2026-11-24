
import sys
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *

Z6 = (Fr(0),)*6
def v6(j, val=1):
    v = list(Z6); v[j] = Fr(val); return tuple(v)
X, Y, P, Q, R, M = (v6(j) for j in range(6))

for name, sym_fn, expect in [
    ("Li31(x,-1) == M P X X", li31_sym, (M, P, X, X)),
    ("Li22(x,-1) == M X P X", li22_sym, (M, X, P, X)),
]:
    S = sym_fn(X, Y, P, Q, R)
    S1 = {}
    for w, cf in S.items():
        w2 = tuple(M if l == R else l for l in w)
        if Y in w2 or Q in w2:
            continue
        S1[w2] = padd(S1.get(w2, {}), cf)
    S1 = {w: cf for w, cf in S1.items() if cf}
    print("PASS" if S1 == {expect: ONE} else "FAIL", name, "" if S1 == {expect: ONE} else S1)
open('progress.log','a').write('step_01: engine validated (y=1,-1 reductions, stuffle, Li4/shuffle killed by d22, d22(F)=COR4 value)\n')
