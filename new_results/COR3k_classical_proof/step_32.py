import sympy as sp, ast
from cor3_common import *

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
CD = load_classical_dict(G1 + "classical_dict.json")
C2 = {
    (0, 1):  sp.polylog(2, x),
    (0, -1): -sp.polylog(2, -x),
    (1, -1): -sp.polylog(2, (1+x)/2) - l1m*l1p + L2*l1p - L2**2/2 + P2/12,
    (-1, 1):  sp.polylog(2, (1+x)/2) - L2*l1p + L2**2/2 - P2/12,
}

def test_word(w_str, nser=3, verbose=True):
    a1, a2, a3 = ast.literal_eval(w_str)
    Cw = CD[w_str]
    Ew = sp.expand(dpoly(Cw) - f_a(a1)*C2[(a2, a3)])
    dE = sp.expand(dpoly(Ew))
    ok1, terms = check_weight1_zero(dE)
    if verbose:
        print(f"(i) dE/dx: {len(terms)} monomials, all rational coefficients 0: {ok1}")
        if not ok1:
            for m, c in terms:
                if c != 0: print("    NONZERO:", m, c)
    c0, resid = const_term(Ew, nser)
    zok, poly = is_zero_const(c0)
    if verbose:
        print(f"(ii) E_w(0+) = {c0}  -> zero polynomial: {zok}")
        print(f"     series residual (must have only x^k, k>=1, times LX^j): {resid}")
    return ok1, zok, resid, Ew

ok1, zok, resid, Ew = test_word('[0, 1, -1]')
