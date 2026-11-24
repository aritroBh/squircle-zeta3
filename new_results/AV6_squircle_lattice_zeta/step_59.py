code = r'''
from Uq import *
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, fsum, mpf
setdps(60)
def zetaEZ_direct_full(s, N, K=1500):
    # sum_{k>j} k^{-(s+N)} j^N directly + tail
    tot = mpf(0)
    for k in range(2, K):
        tot += power(k,-(s+N))*fsum([power(j,N) for j in range(1,k)])
    tail = (zeta(s-1) - fsum([power(k,1-s) for k in range(1,K)]))/(N+1) + (zeta(s)-fsum([power(k,-s) for k in range(1,K)]))/2
    return tot + tail
print("zetaEZ(4,20) formula:", nstr(zetaEZ(4, 20), 25))
print("zetaEZ(4,20) direct :", nstr(zetaEZ_direct_full(4, 20), 25))
print("zetaEZ(4,10) formula:", nstr(zetaEZ(4, 10), 25))
print("zetaEZ(4,10) direct :", nstr(zetaEZ_direct_full(4, 10), 25))
# U1(20) via formula with more n terms
def U1_formula(p, s, nmax):
    tot = mpf(0)
    for n in range(1, nmax+1):
        tot += ((-1)**(n+1))/n * zetaEZ(s, p*n)
    return tot
print("U1(20) formula n=4:", nstr(U1_formula(20, 4, 4), 25))
print("U1(20) formula n=8:", nstr(U1_formula(20, 4, 8), 25))
def U1_direct(p, K=2000):
    s = 4
    tot = mpf(0)
    for k in range(2, K+1):
        inner = mpf(0); kf = mpf(k)
        for j in range(1, k):
            inner += log(1+power(j/kf, p))
        tot += power(kf, -s)*inner
    I0q = fsum([(-1)**(n+1)/(n*(n*p+1)) for n in range(1, 400)])
    from squircle import harmtail
    from mpmath import log as mlog
    tail = I0q*harmtail(K, s-1) - mlog(2)/2*harmtail(K, s)
    return tot + tail
print("U1(20) direct K=2000:", nstr(U1_direct(20, 2000), 25))
print("U1(20) direct K=500 :", nstr(U1_direct(20, 500), 25))
'''
exec(code)
