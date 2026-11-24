import cor8r_sym as S
import sympy as sp
from fractions import Fraction
cc = S.cc
half = sp.Rational(1,2)

d22_cache = {}
def get_Li22(y, conv):
    key=('22',str(y),conv)
    if key not in d22_cache: d22_cache[key]=S.Li22_sym(half, y, conv)
    return d22_cache[key]
def get_Li31(y, conv):
    key=('31',str(y),conv)
    if key not in d22_cache: d22_cache[key]=S.Li31_sym(half, y, conv)
    return d22_cache[key]

# prefactors
p_m = S.LP_add(S.LP_pow_c(-1), S.LP_neg(S.ONE))   # (1-c)/c = c^-1 - 1
p_p = S.LP_add(S.LP_pow_c(-1), S.ONE)             # (1+c)/c = c^-1 + 1

for conv in ['A','B']:
    s22m = get_Li22(1-cc, conv)
    s22p = get_Li22(1+cc, conv)
    s31p = get_Li31(1+cc, conv)
    print(f"conv {conv}: |S Li22(1/2,1-c)|={len(s22m)} terms, |S Li22(1/2,1+c)|={len(s22p)}, |S Li31(1/2,1+c)|={len(s31p)}")
    DT_A = S.S_add(S.S_scaleLP(s22m, p_m), S.S_scaleLP(s22p, S.LP_neg(p_p)))
    DT_B = S.S_add(S.S_scaleLP(s22m, p_m), S.S_scaleLP(s31p, p_p))
    dA = S.delta22(DT_A); dB = S.delta22(DT_B)
    for k in sorted(set(dA)|set(dB)):
        print("   comp", k, " DT_A:", dA.get(k), " DT_B:", dB.get(k))
