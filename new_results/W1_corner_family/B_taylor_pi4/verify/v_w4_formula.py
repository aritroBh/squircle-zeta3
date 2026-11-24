"""v_w4_formula.py -- symbolic (all n) check of the solver's single-sum formula for the Taylor coefficients
w_n of W_4(c) (result.md section 4, typed by the verifier), against the verifier's own canonical Z-sum form
(v_w4_canonical.pkl from v_nested_proof.py).  A = Z_{-1}, B = Z_{-2}, C = Z_{-3}, P = Z_2, Q = Z_3 at 2n;
products are linearised with the verifier's quasi-shuffle.  Also: w_0, w_1, w_2 from direct HPL series.
"""
import pickle, sys
from fractions import Fraction as Fr
import flint
sys.path.insert(0, '.')
from v_nested_proof import RF, X, ONE, emul, acc, load, hpl_series

def E(w):
    return {w: ONE}
def lin(*pairs):
    out = {}
    for cf, Ex in pairs:
        for w, c in Ex.items():
            acc(out, w, (cf if isinstance(cf, RF) else RF(cf))*c)
    return out
def r(num, den):
    return RF(num)/RF(den)

A, B, C, P, Q = E((-1,)), E((-2,)), E((-3,)), E((2,)), E((3,))
N = X
d = 4*N*N - 1
one = {(): ONE}
AB = emul(A, B); AP = emul(A, P); AA = emul(A, A); AAA = emul(AA, A)
rat = lin((r(-(7*N - 4), 4*N**4*(2*N - 1)), one),
          (r(-1, N*d), C),
          (r(8*N**3 + 6*N**2 - 5*N - 2, 2*N**2*(2*N - 1)*(2*N + 1)**2), B),
          (r(3*N - 2, 2*N**3*(2*N - 1)), A),
          (r(6*N + 1, 2*N*(2*N - 1)*(2*N + 1)**2), P),
          (r(2, d), Q),
          (r(1, N*d), lin((1, AB), (-1, Q))),
          (r(-(N - 1), 2*N**2*(2*N - 1)), lin((1, AA), (-1, P))),
          (r(-2, d), lin((1, AP), (-1, C))),
          (r(-1, 6*N*(2*N - 1)), lin((1, AAA), (-3, AP), (2, C))))
z2 = lin((r(-(6*N + 1), 4*N*(2*N - 1)*(2*N + 1)**2), one), (r(1, d), A))
z3 = lin((r(-1, d), one))
# canonical form of [c^(2n)] W_4 recomputed directly with the verifier's nested-sum code
from v_nested_proof import side_hpl
R3 = load('hpl_fit_k3_(0,0,0).json'); PSI = load('hpl_fit_k3_(1,0,0).json'); MPHI = load('hpl_fit_k3_(0,1,0).json')
can = {'1': side_hpl([t for t in R3 if len(t[1]) == 4], 0),
       'z2': side_hpl([t for t in PSI if len(t[1]) == 2], 0),
       'z3': side_hpl([t for t in MPHI if len(t[1]) == 1], 0)}
for part, mine in (('1', rat), ('z2', z2), ('z3', z3)):
    dif = dict(mine)
    for w, c in can[part].items():
        acc(dif, w, -c)
    print(f'w_n single-sum formula ({part}-part) == canonical form of [c^(2n)]W_4, identically in n >= 1: {not dif}')
# low coefficients directly
for nn_ in (0, 1, 2):
    vals = {}
    for part, terms, wt in (('1', R3, 4), ('z2', PSI, 2), ('z3', MPHI, 1)):
        tot = Fr(0)
        for e, w, qv in terms:
            if len(w) == wt and 2*nn_ - e >= 1:
                tot += qv*hpl_series(w, 2*nn_ - e)[2*nn_ - e]
        vals[part] = tot
    print(f'w_{nn_} = {vals["1"]} + ({vals["z2"]}) z2 + ({vals["z3"]}) z3')
