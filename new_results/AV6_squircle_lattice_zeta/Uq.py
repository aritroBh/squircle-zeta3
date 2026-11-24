
"""Exact U_q(p) for integer p via double-zeta + Faulhaber.
U_q(p) = sum_{k>j>=1} k^{-s} ln^q(1+(j/k)^p)
       = q! sum_{n>=q} (-1)^{n-q} |s(n,q)|/n! zeta_EZ(s+pn, -pn),
zeta_EZ(s+N,-N) = (1/(N+1)) sum_{r=0}^N binom(N+1,r) B_r^{(+)} sum_{i=0}^{N+1-r} binom(N+1-r,i)(-1)^i zeta(s-1+r+i).
"""
from mpmath import mp, mpf, zeta, bernoulli, fac, binomial, fsum, log, power, gamma, nstr, fabs

def Bplus(r):
    if r == 1: return mpf(1)/2
    return bernoulli(r)

def zetaEZ(s, N):
    """integer N>=0, s real >=3: sum_{k>j>=1} k^{-(s+N)} j^N."""
    N = int(N)
    tot = mpf(0)
    for r in range(0, N+1):
        inner = mpf(0)
        for i in range(0, N+2-r):
            inner += binomial(N+1-r, i)*((-1)**i)*zeta(s-1+r+i)
        tot += binomial(N+1, r)*Bplus(r)*inner
    return tot/(N+1)

def stirling1u(n, q):
    # unsigned Stirling first kind |s(n,q)|
    if n == q: return 1
    S = [[0]*(q+1) for _ in range(n+1)]
    S[0][0] = 1
    for a in range(1, n+1):
        for b in range(1, min(a,q)+1):
            S[a][b] = S[a-1][b-1] + (a-1)*S[a-1][b]
    return S[n][q]

def Uq_exact(p, s, q, nmax=None):
    """integer p, real s>=3, q>=1."""
    p = int(p); s = mpf(s); q = int(q)
    if nmax is None: nmax = max(q+2, int(20/p)+q+2)
    tot = mpf(0)
    for n in range(q, nmax+1):
        c = ((-1)**(n-q)) * stirling1u(n, q)/fac(n)
        term = c * zetaEZ(s, p*n)
        tot += term
    return fac(q)*tot
