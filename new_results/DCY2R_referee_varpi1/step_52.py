import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, ctx
ctx.dps = 40
from dcy2_cell import *

def gamT3(N, gam):
    D = Delta_of(N); c2m = c2min_of(N)
    T = min(D, c2m/(2*gam))
    return gam*T**3

for (lo, hi) in [(400.0, 408.0), (1000.0, 1020.0), (3075.23, 3136.74), (50000.0, 51000.0), (99000.0, 100000.0)]:
    glo_b3 = B3_of(lo); glo_al = alpha_of(lo)
    ok3 = okA = True
    minr3 = arb(9); minrA = arb(9)
    for k in range(21):
        N = lo + (hi-lo)*k/20
        r3 = gamT3(N, B3_of(N)); rA = gamT3(N, alpha_of(N))
        # cell uses gamma at Nlo and T at Nlo: gamma_lo * T_lo^3 >= gamma(N) T(N)^3 ?
        q3 = gamT3(lo, glo_b3)/r3; qA = gamT3(lo, glo_al)/rA
        if q3 < minr3: minr3 = q3
        if qA < minrA: minrA = qA
    print(f"cell [{lo:9.1f},{hi:9.1f}]: min over cell of cellbound/true  B3: {float(minr3):.4f}  alpha: {float(minrA):.4f}  (>=1 needed)")
