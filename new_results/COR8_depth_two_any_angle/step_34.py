
import sys, json
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
Z3 = (Fr(0),)*3
def v3(j, val=1):
    v = list(Z3); v[j] = Fr(val); return tuple(v)
cv, mv, pv = v3(0), v3(1), v3(2)
zero = v3(0,0)

def drop_zero(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}

S22m = drop_zero(li22_sym(zero, mv, zero, cv, pv))
S22p = drop_zero(li22_sym(zero, pv, zero, cv, mv))
nm = {cv:'c', mv:'m', pv:'p'}
def rd(S):
    return {''.join(nm[l] for l in w): {e: str(v) for e,v in cf.items()} for w,cf in S.items()}
print("S(Li22(1/2,1-c)):", rd(S22m))
print("S(Li22(1/2,1+c)):", rd(S22p))
# also Li31 variants
S31m = drop_zero(li31_sym(zero, mv, zero, cv, pv))
S31p = drop_zero(li31_sym(zero, pv, zero, cv, mv))
print("S(Li31(1/2,1-c)):", rd(S31m))
print("S(Li31(1/2,1+c)):", rd(S31p))
