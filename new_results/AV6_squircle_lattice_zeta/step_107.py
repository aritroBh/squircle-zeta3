code = r'''
from squircle import *
from c2 import G_2, G_l
from mpmath import mp, zeta, mpf, nstr, log, power, fsum, diff as mpdiff
setdps(60)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
s = 4
K = 1500
U1 = mpf(0)
for k in range(1, K+1):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_
        U1 += 4*power(k,-s-1)*L*log(k)
# axes: 2 * k^{-s-1} * (k ln k)(ln k) = 2 k^{-s} ln^2 k
U1 += 2*(zetapp(s) - fsum([power(m,-s)*log(m)**2 for m in range(1,K+1)]))
print("U_1 direct (4x multiplicity):", nstr(U1, 20))
print("8 G_l(4)+4 z''(4):", nstr(8*G_l(4)+4*zetapp(4), 20))
# alternative: 2 x multiplicity version
U1b = mpf(0)
for k in range(1, K+1):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_
        U1b += 2*power(k,-s-1)*L*log(k)
U1b += 2*(zetapp(s) - fsum([power(m,-s)*log(m)**2 for m in range(1,K+1)]))
print("U_1 direct (2x multiplicity):", nstr(U1b, 20))
'''
exec(code)
