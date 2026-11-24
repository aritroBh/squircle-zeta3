import sympy as sp, ast, time
from cor3_common import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
CD = load_classical_dict(G1 + "classical_dict.json")
C2 = {
    (0, 1):  sp.polylog(2, x),
    (0, -1): -sp.polylog(2, -x),
    (1, -1): -sp.polylog(2, (1+x)/2) - l1m*l1p + L2*l1p - L2**2/2 + P2/12,
    (-1, 1):  sp.polylog(2, (1+x)/2) - L2*l1p + L2**2/2 - P2/12,
}

for w_str in sorted(CD.keys()):
    t0 = time.time()
    a1, a2, a3 = ast.literal_eval(w_str)
    Cw = CD[w_str]
    # (i) derivative identity E_w == 0
    Ew = sp.expand(dpoly(Cw) - f_a(a1)*C2[(a2, a3)])
    okE, termsE = check_weight2_zero(Ew)
    # (ii) lim_{x->0+} C_w = 0
    c0, resid = const_term(Cw, 2)
    okc, poly = is_zero_const(c0)
    status = "PROVED" if (okE and okc) else "FAILED"
    print(f"w = {w_str:>12}: derivative identity EXACT: {okE}; constant EXACT 0: {okc}"
          f" (c0 = {c0}) -> {status}   [{time.time()-t0:.1f}s]")
    if not okE:
        for m, c in termsE:
            if c != 0: print("      NONZERO:", m, "=", c)
