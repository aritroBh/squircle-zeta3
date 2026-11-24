"""v_flat_numeric.py -- INDEPENDENT numerical check of Theorem A (verifier's code).

Flat-angle expansion: M(u;alpha) = -sum_n C(2n,n)/(2n-1) c^(2n) J_n(u) with
   J_n(u) = int_0^1 (1+g)^(-1-u) h^n (1+h)^(1-2n) g^(u-1) dg,  h = g^(1-u)
(ORIGINAL g-variable on [0,1]; no use of the g -> 1/g symmetry or the y-substitution or Beta derivatives).
Taylor coefficients in u are computed by truncated power-series arithmetic at each quadrature node.
mu_{k,n} = -C(2n,n)/(2n-1) [u^k] J_n.
Checks: (1) agreement with the solver's exact table flat_exact_N40_K3.json;
        (2) PSLQ of mu_{2,n} on {1, z2, log2, log2^2} and mu_{3,n} on {1, z2, z3, log2, log2^2, log2^3, z2 log2}.
"""
import json
import mpmath as mp
from fractions import Fraction as Fr

mp.mp.dps = 70
K = 3


def smul(A, B):
    return [mp.fsum(A[i]*B[m - i] for i in range(m + 1)) for m in range(K + 1)]


def sexp(A):          # exp of series
    E = [mp.exp(A[0])] + [mp.mpf(0)]*K
    for m in range(1, K + 1):
        E[m] = mp.fsum(i*A[i]*E[m - i] for i in range(1, m + 1))/m
    return E


def slog(A):          # log of series with A[0] > 0
    Lg = [mp.log(A[0])] + [mp.mpf(0)]*K
    for m in range(1, K + 1):
        Lg[m] = (A[m] - mp.fsum(i*Lg[i]*A[m - i] for i in range(1, m))/m)/A[0]
    return Lg


def integrand(g, n):
    L = mp.log(g); m1 = mp.log1p(g)
    # h = g e^{-uL}
    h = [g*(-L)**k/mp.factorial(k) for k in range(K + 1)]
    logh = [L, -L] + [mp.mpf(0)]*(K - 1)
    one_h = [1 + h[0]] + h[1:]
    sl = slog(one_h)
    lg = [n*logh[k] + (1 - 2*n)*sl[k] + (-(1 + 0)*m1 if k == 0 else 0) for k in range(K + 1)]
    # (1+g)^(-1-u): log = -(1+u) m1 ; g^(u-1): log = (u-1)L
    lg[0] += -L
    lg[1] += -m1 + L
    return sexp(lg)


def J(n):
    cache = {}

    def F(g, k):
        if g not in cache:
            cache[g] = integrand(g, n)
        return cache[g][k]
    return [mp.quad(lambda g: F(g, k), [0, mp.mpf(1)/16, mp.mpf(1)/2, 1]) for k in range(K + 1)]


tab = json.load(open('flat_exact_N40_K3.json'))
z2, z3 = mp.zeta(2), mp.zeta(3)
L2 = mp.log(2)
maxdiff = 0
for n in range(1, 9):
    Jn = J(n)
    pref = mp.binomial(2*n, n)/(2*n - 1)
    mu = [-pref*x for x in Jn]
    for k in range(K + 1):
        e = tab[str(k)][str(n)]
        val = mp.mpf(0)
        for key, s in e.items():
            fr = Fr(s); mon = eval(key)
            val += mp.mpf(fr.numerator)/fr.denominator*z2**mon[0]*z3**mon[1]
        d = abs(val - mu[k])
        maxdiff = max(maxdiff, d)
        line = f'n={n} k={k} mu={mp.nstr(mu[k], 20)} |mu - exact table| = {mp.nstr(d, 3)}'
        if k == 2:
            rel = mp.pslq([mu[k], 1, z2, L2, L2**2], maxcoeff=10**12, maxsteps=10**6, tol=mp.mpf(10)**-50)
            line += f'  PSLQ[mu,1,z2,log2,log2^2]={rel}'
        if k == 3:
            rel = mp.pslq([mu[k], 1, z2, z3, L2, L2**2, L2**3, z2*L2], maxcoeff=10**14, maxsteps=10**6, tol=mp.mpf(10)**-50)
            line += f'  PSLQ[mu,1,z2,z3,log2,log2^2,log2^3,z2log2]={rel}'
        print(line, flush=True)
print('max |numeric - exact table| over n<=8, k<=3:', mp.nstr(maxdiff, 3))
