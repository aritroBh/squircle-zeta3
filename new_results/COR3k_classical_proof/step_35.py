import sympy as sp
from cor3_common import *

# Classical weight-2 atom identities (to be proved below):
LI2_REDUCE = {
    'x': sp.polylog(2, x),
    '-x': sp.polylog(2, -x),
    '1-x': -sp.polylog(2, x) + P2/6 - logx*l1m,
    '1/(1+x)': sp.polylog(2, -x) + P2/6 + logx*l1p - l1p**2/2,
    '(1+x)/2': sp.polylog(2, (1+x)/2),
    '(1-x)/(1+x)': sp.polylog(2, (1-x)/(1+x)),
    '-(1-x)/(1+x)': sp.polylog(2, (1+x)/2) - P2/6 + l1m*l1p - l1m*L2 - l1p**2/2 + L2**2/2,
    '2x/(1+x)': -sp.polylog(2, (1-x)/(1+x)) + P2/6 - (L2 + logx - l1p)*(l1m - l1p),
}

print("=== PROOF of the 4 non-trivial Li2 atom identities ===\n")
for nm in ['1-x', '1/(1+x)', '-(1-x)/(1+x)', '2x/(1+x)']:
    lhs = sp.polylog(2, ARGS[nm])
    G = sp.expand(lhs - LI2_REDUCE[nm])          # claimed identically 0
    dG = sp.expand(dpoly(G))
    okd, terms = check_weight1_zero(dG)
    c0, resid = const_term(G, 3)
    okc, poly = is_zero_const(c0)
    print(f"Li2({nm}): derivative log-coeffs all 0: {okd};  G(0+) = {c0} zero-poly: {okc}"
          f"  -> {'PROVED' if okd and okc else 'FAILED'}")
    if not okd:
        for m, c in terms:
            if c != 0: print("   NONZERO:", m, c)
