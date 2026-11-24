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

def prove_word_E(w_str, verbose=True):
    a1, a2, a3 = ast.literal_eval(w_str)
    Cw = CD[w_str]
    Ew = sp.expand(dpoly(Cw) - f_a(a1)*C2[(a2, a3)])
    ok, terms = check_weight2_zero(Ew)
    if verbose:
        print(f"  E_w after Li2 reduction: {len(terms)} monomials; all rational coefficients 0: {ok}")
        if not ok:
            for m, c in terms:
                if c != 0: print("    NONZERO:", m, "=", c)
    return ok, Ew

t0 = time.time()
print("=== [0, 1, -1]: prove E_w == 0 (derivative identity) ===")
ok, Ew = prove_word_E('[0, 1, -1]')
print(f"  time {time.time()-t0:.1f}s")
