"""flat_expansion.py -- Taylor expansion of the corner coefficients about the FLAT angle alpha = pi.

Exact expansion (derivation in result.md, Sec. 4):  with h = g^(1-u) and c = cos(alpha/2),
  1 + h^2 - 2 cos(alpha) h = (1+h)^2 - 4 c^2 h,
  (1+h) sqrt(1 - 4c^2 h/(1+h)^2) = (1+h) - sum_{n>=1} kappa_n 4^n c^(2n) h^n (1+h)^(1-2n),
  kappa_n = binom(2n,n)/((2n-1) 4^n)   [sqrt(1-z) = 1 - sum kappa_n z^n],
and M(u; pi) = 0 identically, so
  M(u; alpha) = - sum_{n>=1} kappa_n 4^n c^(2n) J_n(u),
  J_n(u) = int_0^1 (1+g)^(-1-u) g^(n(1-u)+u-1) (1+g^(1-u))^(1-2n) dg   (converges for c^2 < 1).
The coefficient of c^(2n) in M_k is mu_{k,n} = -kappa_n 4^n [u^k] J_n(u).
mu_{k,n} is an alternating MZV combination of weight <= k+1.  We PSLQ each mu_{k,n} against the basis of
alternating MZVs of weight <= k+1 (real, letters 0,-1):
  w0: 1; w1: log2; w2: pi^2, log^2 2; w3: zeta3, pi^2 log2, log^3 2;
  w4: Li4(1/2), zeta3 log2, pi^4, pi^2 log^2 2, log^4 2.
and report which weights occur.  Also checks sum_n mu_{k,n} c^(2n) against stored M_k at 2pi/3 (c^2=1/4).
"""
import sys
import mpmath as mp

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 160
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 8
KMAX = 3
mp.mp.dps = DPS + 20


def smul(A, B, N):
    return [mp.fsum(A[i]*B[n - i] for i in range(n + 1)) for n in range(N + 1)]


def sexp(A, N):          # exp of a series with A[0] arbitrary
    e0 = mp.exp(A[0])
    B = [mp.mpf(1)] + [mp.mpf(0)]*N
    # B' = A' B
    for n in range(1, N + 1):
        B[n] = mp.fsum(k*A[k]*B[n - k] for k in range(1, n + 1))/n
    return [e0*b for b in B]


def slog1p_of(Y, N):     # log(1+Y), Y[0] > -1
    # (1+Y) L' = Y'
    one = [1 + Y[0]] + Y[1:]
    Lg = [mp.log(one[0])] + [mp.mpf(0)]*N
    for n in range(1, N + 1):
        s = n*Y[n] - mp.fsum(k*Lg[k]*one[n - k] for k in range(1, n))
        Lg[n] = s/(n*one[0])
    return Lg


def Jn_series(g, n, N):
    L = mp.log(g)
    m = mp.log1p(g)
    # log of integrand, as a series in u (all factors are exp(...))
    # (1+g)^(-1-u): -(1+u) m ; g^(n(1-u)+u-1): (n-1)L - (n-1)uL ; (1+g e^{-uL})^(1-2n)
    Y = [g*(-L)**k/mp.factorial(k) for k in range(N + 1)]
    l1y = slog1p_of(Y, N)
    lg = [mp.mpf(0)]*(N + 1)
    lg[0] = -m + (n - 1)*L + (1 - 2*n)*l1y[0]
    if N >= 1:
        lg[1] = -m - (n - 1)*L
    for k in range(1, N + 1):
        lg[k] += (1 - 2*n)*l1y[k]
    return sexp(lg, N)


def Jn_coeffs(n, N):
    cache = {}

    def F(g, k):
        if g not in cache:
            cache[g] = Jn_series(g, n, N)
        return cache[g][k]
    pts = [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1]
    return [mp.quad(lambda g: F(g, k), pts) for k in range(N + 1)]


def basis(wmax):
    l2, z3, pi = mp.log(2), mp.zeta(3), mp.pi
    B = [('1', mp.mpf(1)), ('log2', l2), ('pi^2', pi**2), ('log2^2', l2**2), ('zeta3', z3),
         ('pi^2log2', pi**2*l2), ('log2^3', l2**3), ('Li4(1/2)', mp.polylog(4, mp.mpf(1)/2)),
         ('zeta3log2', z3*l2), ('pi^4', pi**4), ('pi^2log2^2', pi**2*l2**2), ('log2^4', l2**4)]
    nb = {0: 1, 1: 2, 2: 4, 3: 7, 4: 12}[wmax]
    return B[:nb]


if __name__ == '__main__':
    from fractions import Fraction
    rows = {}
    for n in range(1, NMAX + 1):
        kap = mp.binomial(2*n, n)/((2*n - 1)*mp.mpf(4)**n)
        J = Jn_coeffs(n, KMAX)
        rows[n] = [-kap*4**n*J[k] for k in range(KMAX + 1)]
        for k in range(KMAX + 1):
            B = basis(min(k + 1, 4))
            vec = [rows[n][k]] + [b[1] for b in B]
            mp.mp.dps = DPS
            rel = mp.pslq(vec, maxcoeff=10**30, maxsteps=10**6)
            mp.mp.dps = DPS + 20
            if rel is None:
                print(f'n={n} k={k}: mu = {mp.nstr(rows[n][k], 20)}  NO RELATION', flush=True)
                continue
            d = rel[0]
            terms = {B[i][0]: Fraction(-rel[i + 1], d) for i in range(len(B)) if rel[i + 1] != 0}
            resid = abs(mp.fsum(r*v for r, v in zip(rel, vec)))
            print(f'n={n} k={k}: mu = {mp.nstr(rows[n][k], 20)} = {terms}  '
                  f'(maxcoef {max(abs(r) for r in rel)}, resid {mp.nstr(resid, 3)})', flush=True)
    # partial-sum sanity check at 2pi/3 (c^2 = 1/4) against stored values
    vals = {}
    for line in open('values_2pi_3.txt'):
        if line.startswith('M'):
            k, v = line.split('=')
            vals[int(k.strip()[1:])] = mp.mpf(v.strip())
    for k in range(KMAX + 1):
        ps = mp.fsum(rows[n][k]*mp.mpf(1)/4**n for n in rows)
        print(f'2pi/3: M{k} - partial sum (n<={NMAX}) = {mp.nstr(vals[k] - ps, 5)}  (next term ~ {mp.nstr(rows[NMAX][k]/4**NMAX, 3)})')
