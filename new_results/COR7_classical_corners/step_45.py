import math
import mpmath as mp

def _mul(A, B, N):
    return [mp.fsum(A[i]*B[n - i] for i in range(n + 1)) for n in range(N + 1)]

def _sqrt(Q, N):
    b = [mp.sqrt(Q[0])]
    for n in range(1, N + 1):
        s = Q[n] - mp.fsum(b[k]*b[n - k] for k in range(1, n))
        b.append(s/(2*b[0]))
    return b

def integrand_series_hp(g, c, N):
    """u-series of the corner integrand at node g, with internal precision boosted
    by the proximity of g to 0 (fixes the (1+g)^-1*sqrt(...)-1 cancellation)."""
    base = mp.mp.dps
    boost = 0
    if g < mp.mpf('0.1'):
        boost = int(-mp.log10(g)) + 15
    mp.mp.dps = base + boost
    g2 = mp.mpf(g)  # renormalize to current precision
    c2 = mp.mpf(c)
    L = mp.log(g2)
    m = mp.log1p(g2)
    fac = [mp.factorial(k) for k in range(N + 1)]
    E1 = [(-L)**k/fac[k] for k in range(N + 1)]
    E2 = [(-2*L)**k/fac[k] for k in range(N + 1)]
    Q = [1 + g2*g2 - 2*c2*g2] + [g2*g2*E2[k] - 2*c2*g2*E1[k] for k in range(1, N + 1)]
    B = _sqrt(Q, N)
    A = [(-m)**k/fac[k]/(1 + g2) for k in range(N + 1)]
    AB = _mul(A, B, N)
    AB[0] -= 1
    G = [L**k/fac[k]/g2 for k in range(N + 1)]
    out = _mul(AB, G, N)
    mp.mp.dps = base
    return out

def corner_coeffs_hp(a, N, pts=None):
    c = mp.cos(a)
    cache = {}
    def F(g, k):
        key = g
        if key not in cache:
            cache[key] = integrand_series_hp(g, c, N)
        return cache[key][k]
    if pts is None:
        pts = [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1]
    return [mp.quad(lambda g: F(g, k), pts) for k in range(N + 1)]

import sys
sys.path.insert(0, '.')
import cor7_lib as L
# test at c=1/3 against HPL formula
for dps in (100, 200, 300):
    mp.mp.dps = dps
    c3 = mp.mpf(1)/3
    M3ref, W4, W3, W2 = L.M3_weights(c3)
    M3q = corner_coeffs_hp(2*mp.acos(c3), 3)[3]
    print(f"dps {dps}: |quad_hp - HPL| = {mp.nstr(abs(M3q - M3ref), 3)}")
