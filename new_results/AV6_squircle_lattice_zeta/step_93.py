code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum, diff as mpdiff
import mpmath
setdps(60)
s = 4
def zetap(a): return mpdiff(lambda x: zeta(x), a)
def HTln(sig, K):
    return -zetap(sig) - fsum([power(m,-sig)*log(m) for m in range(1, K+1)])
def HT(sig, K):
    return zeta(sig) - fsum([power(m,-sig) for m in range(1, K+1)])
K = 2000
A = mpf(0); B = mpf(0)
for k in range(1, K+1):
    sig = fsum([m0*log(m0) for m0 in range(1, k)])
    A += 4*k*power(k,-s)*log(k)
    B += power(k,-s-1)*(4*sig + 2*k*log(k))
A += 4*HTln(s-1, K)
lnA = log(mpmath.glaisher)
Bt = mpf(0)
Bt += 4*( mpf(1)/2*HTln(s-1, K) - mpf(1)/4*HT(s-1, K) - mpf(1)/2*HTln(s, K) + mpf(1)/12*HTln(s+1, K) + lnA*HT(s+1, K) )
for i in range(1, 40):
    Bt += 4*( -mpf(1)/(2*i)*HT(s-1+i, K) + mpf(1)/(2*i)*HT(s+i, K) - mpf(1)/(12*i)*HT(s+1+i, K) )
Bt += 2*HTln(s, K)
B += Bt
A *= s; B *= s
print("A =", nstr(A, 25), " (expect -4s zeta'(3) =", nstr(-4*s*zetap(s-1), 15), ")")
print("B =", nstr(B, 25))
print("c1 = A - B =", nstr(A-B, 25))
print("FD value: ~1.247")
'''
exec(code)
