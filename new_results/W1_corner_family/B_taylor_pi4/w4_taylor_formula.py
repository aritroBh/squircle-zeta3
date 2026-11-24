"""w4_taylor_formula.py -- PROVED (exact): the flat-angle Taylor coefficients of the weight-4 component
W_4(c) = sum_{n>=0} w_n c^{2n} of M_3(alpha), written with single sums only, and checked exactly against
(i) the canonical Z-sum form from b3_allN_proof.py (w4_taylor_zsum.pkl), and (ii) the exact Taylor coefficients
of the weight-4 part of B3 (HPL series in exact rationals), for 1 <= n <= NCHK.

Notation (all at argument 2n):
  A = Z_{-1} = H_n - H_{2n},  B = Z_{-2} = H^(2)_n/2 - H^(2)_{2n},  C = Z_{-3} = H^(3)_n/4 - H^(3)_{2n},
  P = Z_2 = H^(2)_{2n},  Q = Z_3 = H^(3)_{2n}.
Uses the elementary-symmetric identities (x_i = (-1)^i/i):  Z_{-1,-1} = (A^2 - P)/2,
Z_{-1,-1,-1} = (A^3 - 3AP + 2C)/6,  Z_{-2,-1} + Z_{-1,-2} = AB - Q,  Z_{-1,2} + Z_{2,-1} = AP - C.
"""
import sys, pickle
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, '.')
sys.setrecursionlimit(10000)
import zsum
from hpl import coeffs_exact
from check_b3_split import R3, PSI, MPHI

NCHK = int(sys.argv[1]) if len(sys.argv) > 1 else 150
n = sp.Symbol('n')


def H(r, m):
    return sum((Fr(1, i**r) for i in range(1, m + 1)), Fr(0))


def single_form(nv):
    A = H(1, nv) - H(1, 2*nv); B = H(2, nv)/2 - H(2, 2*nv); C = H(3, nv)/4 - H(3, 2*nv)
    P = H(2, 2*nv); Q = H(3, 2*nv)
    N = Fr(nv)
    d = 4*N*N - 1
    rat = (-(7*N - 4)/(4*N**4*(2*N - 1))
           - C/(N*d)
           + (8*N**3 + 6*N**2 - 5*N - 2)/(2*N**2*(2*N - 1)*(2*N + 1)**2)*B
           + (3*N - 2)/(2*N**3*(2*N - 1))*A
           + (6*N + 1)/(2*N*(2*N - 1)*(2*N + 1)**2)*P
           + 2*Q/d
           + (A*B - Q)/(N*d)
           - (N - 1)/(N**2*(2*N - 1))*(A*A - P)/2
           - 2*(A*P - C)/d
           - (A**3 - 3*A*P + 2*C)/(6*N*(2*N - 1)))
    z2 = -(6*N + 1)/(4*N*(2*N - 1)*(2*N + 1)**2) + A/d
    z3 = -1/d
    return {'1': rat, 'z2': z2, 'z3': z3}


def from_zsum(nv, E):
    tot = Fr(0)
    for w, c in E:
        cv = sp.Rational(sp.sympify(c).subs(n, nv))
        tot += Fr(int(cv.p), int(cv.q))*zsum.Z_val(tuple(w), 2*nv)
    return tot


def hpl_part(terms, wt, Nmax):
    """exact Taylor coefficients (index = power of c) of sum v c^e H_w over |w| = wt."""
    out = [Fr(0)]*(Nmax + 1)
    for e, w, v in terms:
        if len(w) != wt:
            continue
        ce = coeffs_exact(w, Nmax + 3)
        for N in range(0, Nmax + 1):
            idx = N - e
            if 1 <= idx <= Nmax + 2:
                out[N] += v*ce[idx]
    return out


if __name__ == '__main__':
    Z = pickle.load(open('w4_taylor_zsum.pkl', 'rb'))
    Nmax = 2*NCHK + 1
    hp = {'1': hpl_part(R3, 4, Nmax), 'z2': hpl_part(PSI, 2, Nmax), 'z3': hpl_part(MPHI, 1, Nmax)}
    ok1 = ok2 = True
    for nv in range(1, NCHK + 1):
        sf = single_form(nv)
        for part in ('1', 'z2', 'z3'):
            if nv <= 60 and sf[part] != from_zsum(nv, Z[part]):
                ok1 = False; print('single-sum form != Z-sum form', part, nv)
            if sf[part] != hp[part][2*nv]:
                ok2 = False; print('single-sum form != HPL Taylor coefficient', part, nv)
    odd_zero = all(hp[p][2*k + 1] == 0 for p in hp for k in range(0, NCHK))
    print(f'single-sum form == Z-sum canonical form (n<=60): {ok1}')
    print(f'single-sum form == exact Taylor coefficients of the weight-4 part of B3 (n<={NCHK}): {ok2}')
    print('odd coefficients of the weight-4 part vanish (exact, c^1..c^%d):' % (2*NCHK - 1), odd_zero)
    print('constant term w_0 = W_4(0): rational part', hp['1'][0], '; z2 part', hp['z2'][0], '; z3 part', hp['z3'][0])
    print('w_1 =', single_form(1), '  w_2 =', single_form(2))
