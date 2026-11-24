code = r'''
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
'''
with open("Uq.py","w") as f: f.write(code)
test = r'''
from Uq import *
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, fsum, mpf
setdps(60)
# check zetaEZ against direct sum
def zetaEZ_direct(s, N, K=3000):
    tot = mpf(0)
    for k in range(2, K):
        tot += power(k,-(s+N))*fsum([power(j,N) for j in range(1,k)])
    return tot
s = 4
print("zetaEZ(5,1):", nstr(zetaEZ(5,1), 25))
print("direct     :", nstr(zetaEZ_direct(5,1), 25))
print("zetaEZ(6,2):", nstr(zetaEZ(6,2), 25))
print("direct     :", nstr(zetaEZ_direct(6,2), 25))
# check Uq_exact against U1_direct
def U1_direct(p, K=3000):
    tot = mpf(0)
    for k in range(2, K+1):
        inner = mpf(0)
        kf = mpf(k)
        for j in range(1, k):
            inner += log(1+power(j/kf, p))
        tot += power(kf, -s)*inner
    I0q = fsum([(-1)**(n+1)/(n*(n*p+1)) for n in range(1, 400)])
    from squircle import harmtail
    from mpmath import log as mlog
    tail = I0q*harmtail(K, s-1) - mlog(2)/2*harmtail(K, s)
    return tot + tail
for p_ in [5, 10, 20]:
    print(f"U1({p_}): exact={nstr(Uq_exact(p_, s, 1), 25)}")
    print(f"        direct={nstr(U1_direct(p_), 25)}")
'''
exec(test)
with open("progress.log","a") as f: f.write("step_40: Uq.py exact U_q via double zeta\n")
