code = r'''
from Uq import *
from squircle import setdps, harmtail
from mpmath import mp, nstr, fabs, log, power, fsum, mpf
setdps(60)
def zetaEZ_direct(s, N, K=2000):
    tot = mpf(0)
    for k in range(2, K):
        tot += power(k,-(s+N))*fsum([power(j,N) for j in range(1,k)])
    # tail: H_{k-1}^{(-N)} ~ k^{N+1}/(N+1)+k^N/2: tail ~ zeta(s-1)-H...
    tail = (zeta(s-1) - fsum([power(k,1-s) for k in range(1,K)]))/(N+1) + (zeta(s)-fsum([power(k,-s) for k in range(1,K)]))/2
    return tot + tail
print("zetaEZ(24,-20) exact :", nstr(zetaEZ(24,20), 25))
print("zetaEZ(24,-20) direct:", nstr(zetaEZ_direct(24,20), 25))
print("zetaEZ(14,-10) exact :", nstr(zetaEZ(14,10), 25))
print("zetaEZ(14,-10) direct:", nstr(zetaEZ_direct(14,10), 25))
'''
exec(code)
