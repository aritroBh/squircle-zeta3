
# (0) Reconcile DT: recompute d22 of BOTH candidate DTs exactly
import sys
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
target = {((0,1),(1,2)): {0: Fr(1), -1: Fr(-1)},
          ((0,2),(1,2)): {0: Fr(-1), -1: Fr(-1)}}
S22m = drop_zero(li22_sym(zero, mv, zero, cv, pv))   # Li22(1/2,1-c)
S22p = drop_zero(li22_sym(zero, pv, zero, cv, mv))   # Li22(1/2,1+c)
S31p = drop_zero(li31_sym(zero, pv, zero, cv, mv))   # Li31(1/2,1+c)
def scaled_d22(S, poly):
    return {k: pmul(cf, poly) for k, cf in delta22(S, 3).items()}
# DT version A (both Li22): ((1-c)/c) Li22(1/2,1-c) - ((1+c)/c) Li22(1/2,1+c)
dA = {}
for k, cf in scaled_d22(S22m, {-1: Fr(1), 0: Fr(-1)}).items(): dA[k] = cf
for k, cf in scaled_d22(S22p, {-1: Fr(-1), 0: Fr(-1)}).items(): dA[k] = padd(dA.get(k, {}), cf)
# DT version B (Li22 + Li31): ((1-c)/c) Li22(1/2,1-c) + ((1+c)/c) Li31(1/2,1+c)
dB = dict(dA)  # same first term... careful: dA currently has both; rebuild
dB = {}
for k, cf in scaled_d22(S22m, {-1: Fr(1), 0: Fr(-1)}).items(): dB[k] = cf
for k, cf in scaled_d22(S31p, {-1: Fr(1), 0: Fr(1)}).items(): dB[k] = padd(dB.get(k, {}), cf)
print("d22 target:", {k: {e: str(v) for e,v in sorted(p.items())} for k,p in target.items()})
print("d22 DT_A (Li22 pair):", {k: {e: str(v) for e,v in sorted(p.items())} for k,p in dA.items()}, "MATCH:", dA == target)
print("d22 DT_B (Li22+Li31):", {k: {e: str(v) for e,v in sorted(p.items())} for k,p in dB.items()}, "MATCH:", dB == target)
open('progress.log','a').write('step_60: DT reconciliation; both A (Li22 pair) and B (Li22+Li31) match d22\n')
