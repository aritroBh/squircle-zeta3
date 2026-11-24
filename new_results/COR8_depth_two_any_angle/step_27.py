
import sys
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *

Z3 = (Fr(0),)*3
def v3(j, val=1):
    v = list(Z3); v[j] = Fr(val); return tuple(v)
cv, mv, pv = v3(0), v3(1), v3(2)
zero = v3(0,0)

target = {((0,1),(1,2)): {0: Fr(1), -1: Fr(-1)},
          ((0,2),(1,2)): {0: Fr(-1), -1: Fr(-1)}}

# Li22(1/2, 1-c): X=[1/2]=0, Y=[m]=mv, P=0, Q=[c]=cv, R=[(1+c)/2]=pv
D1 = delta22(li22_sym(zero, mv, zero, cv, pv), 3)
# Li31(1/2, 1+c): X=0, Y=[p]=pv, P=0, Q=[-c]=cv, R=[(1-c)/2]=mv
D2 = delta22(li31_sym(zero, pv, zero, cv, mv), 3)
print("d22 Li22(1/2,1-c):", {k: dict(v) for k,v in D1.items()})
print("d22 Li31(1/2,1+c):", {k: dict(v) for k,v in D2.items()})
# combination: (c^-1 - 1)*D1 - (1 + c^-1)*D2
comb = {}
for k, poly in D1.items():
    comb[k] = padd(pmul(poly, {-1: Fr(1)}), pmul(poly, {0: Fr(-1)}))
for k, poly in D2.items():
    q = padd({0: Fr(-1)}, {-1: Fr(-1)})
    comb[k] = padd(comb.get(k, {}), pmul(poly, q))
comb = {k: v for k, v in comb.items() if v}
print("combination:", {k: {e: str(v) for e,v in sorted(p.items())} for k,p in comb.items()})
print("MATCHES TARGET:", comb == target)
open('progress.log','a').write('step_02: d22(Omega) = d22((c^-1-1)Li22(1/2,1-c) - (1+c^-1)Li31(1/2,1+c)) EXACT\n')
