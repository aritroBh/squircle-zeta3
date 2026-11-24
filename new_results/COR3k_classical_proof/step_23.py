import sympy as sp
from cor3_common import *

# Weight-1 HPLs (definitions): H_0 = log x, H_1 = -log(1-x), H_{-1} = log(1+x)
H1w = {0: logx, 1: -l1m, -1: l1p}

# Weight-2 classical forms (proved here):
C2 = {
    (0, 1):  sp.polylog(2, x),
    (0, -1): -sp.polylog(2, -x),
    (1, 1):  l1m**2/2,
    (-1, -1): l1p**2/2,
    (1, -1): -sp.polylog(2, (1+x)/2) - l1m*l1p + L2*l1p - L2**2/2 + P2/12,
    (-1, 1):  sp.polylog(2, (1+x)/2) - L2*l1p + L2**2/2 - P2/12,
}

print("=== WEIGHT 2: prove C_{a,b} = H_{a,b} on (0,1) ===\n")
for (a, b), Cab in C2.items():
    # (i) d/dx C_{a,b} - f_a * H_b  in the log basis
    D = sp.expand(dpoly(Cab) - f_a(a)*H1w[b])
    ok, coeffs = check_weight1_zero(D)
    # (ii) value at 0+
    c0, resid = const_term(Cab, 3)
    zok, poly = is_zero_const(c0)
    print(f"H[{a},{b}]: derivative-log-coeffs all 0: {ok} {coeffs}")
    print(f"        C(0+) = {sp.expand(c0)} -> zero polynomial: {zok}; series residual: {resid}")
    print(f"        -> {'PROVED' if (ok and zok) else 'FAILED'}\n")
