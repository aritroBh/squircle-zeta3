"""
H1: high-order Taylor coefficients (numerical, high precision) of S_n(1+d)/S_n(1), via Route A
(same exact formula as 01_series_routeA.py, but with numerical zeta values so that K can be large).

  S_n/S_n(1) = (n-1)! eps^{n-1}/Gamma(n-(n-1)eta) * sum_{k<=K} binom(1/2,k) n^{-k} T_k,
  T_k = k! [s^k] Phi(s)^n,  Phi(s) = sum_k Phi_k s^k/k!,  Phi_k = sum_j C(k,j)(-1)^{k-j} Gamma(1+(2j-1)eta).
Usage: python 05_high_order_numeric.py n K   -> writes hi_coeffs_n{n}_K{K}.json (eta- and d-coefficients)
"""
import sys, json, math
import mpmath as mp

n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
K = int(sys.argv[2]) if len(sys.argv) > 2 else 40
mp.mp.dps = 80 + int(3.2 * K)


def ser_mul(a, b, lo=0):
    """truncated product of eta-series (lists of length K+1)"""
    c = [mp.mpf(0)] * (K + 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j in range(0, K + 1 - i):
            bj = b[j]
            if bj != 0:
                c[i + j] += ai * bj
    return c


def ser_exp(L):
    """exp of series L with L[0] = 0: E' = L' E"""
    E = [mp.mpf(0)] * (K + 1)
    E[0] = mp.mpf(1)
    dL = [k * L[k] for k in range(K + 1)]
    for k in range(1, K + 1):
        s = mp.mpf(0)
        for j in range(1, k + 1):
            s += dL[j] * E[k - j]
        E[k] = s / k
    return E


Lg = [mp.mpf(0), -mp.euler] + [(-1) ** j * mp.zeta(j) / j for j in range(2, K + 1)]


def gamma1_ser(c):
    """Gamma(1 + c*eta)"""
    return ser_exp([Lg[r] * mp.mpf(c) ** r for r in range(K + 1)])


G = [gamma1_ser(2 * j - 1) for j in range(K + 1)]
Phi = []
for k in range(K + 1):
    s = [mp.mpf(0)] * (K + 1)
    for j in range(k + 1):
        w = math.comb(k, j) * (-1) ** (k - j)
        for r in range(K + 1):
            s[r] += w * G[j][r]
    for r in range(k):  # Phi_k = O(eta^k); these entries are pure cancellation noise
        assert abs(s[r]) < mp.mpf(10) ** (-mp.mp.dps + 3 * K + 40), (k, r, s[r])
        s[r] = mp.mpf(0)
    Phi.append(s)

# EGF power: P = Phi^n
P = [row[:] for row in Phi]
for _ in range(n - 1):
    Q = []
    for k in range(K + 1):
        acc = [mp.mpf(0)] * (K + 1)
        for i in range(k + 1):
            prod = ser_mul(P[i], Phi[k - i])
            w = math.comb(k, i)
            for r in range(k, K + 1):
                acc[r] += w * prod[r]
        Q.append(acc)
    P = Q
E = [mp.mpf(0)] * (K + 1)
for k in range(K + 1):
    bk = mp.binomial(mp.mpf(1) / 2, k) / mp.mpf(n) ** k
    for r in range(K + 1):
        E[r] += bk * P[k][r]
# prefactor (1-eta)^{n-1} / Gamma(1-(n-1)eta) / prod_{i=1}^{n-1}(1-(n-1)eta/i)
pref = ser_exp([-Lg[r] * mp.mpf(-(n - 1)) ** r for r in range(K + 1)])
for i in range(1, n):
    geo = [(mp.mpf(n - 1) / i) ** r for r in range(K + 1)]
    pref = ser_mul(pref, geo)
for i in range(n - 1):
    pref = ser_mul(pref, [mp.mpf(1), mp.mpf(-1)] + [mp.mpf(0)] * (K - 1))
A = ser_mul(pref, E)  # coefficients in eta
# to d: eta = d/(1+d)
etad = [mp.mpf(0)] + [mp.mpf((-1) ** (m - 1)) for m in range(1, K + 1)]
B = [mp.mpf(0)] * (K + 1)
pw = [mp.mpf(1)] + [mp.mpf(0)] * K
for e in range(K + 1):
    if e > 0:
        pw = ser_mul(pw, etad)
    for r in range(K + 1):
        B[r] += A[e] * pw[r]
out = {'n': n, 'K': K, 'eta': [mp.nstr(x, 150) for x in A], 'd': [mp.nstr(x, 150) for x in B]}
json.dump(out, open(f'hi_coeffs_n{n}_K{K}.json', 'w'), indent=1)
mp.mp.dps = 30
for k in range(K + 1):
    print(k, mp.nstr(A[k], 25), mp.nstr(B[k], 25))
print('saved', f'hi_coeffs_n{n}_K{K}.json')
