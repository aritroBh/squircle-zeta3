code = r'''
from squircle import *
from c2 import G_2, G_l, G_x, G_22
from mpmath import mp, zeta, mpf, nstr, log, power, fsum, diff as mpdiff
setdps(60)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
s = 3
K = 1500
# U_1 = 8 G_l(s) + 4 zeta''(s)
U1d = 8*G_l(3) 
print("8 G_l(3):", nstr(U1d, 20), " +4 z''(3) =", nstr(U1d+4*zetapp(3), 20))
# U_2 = 8(G_22+G_x)(s+1)+4 z''(s)
print("8(G_22+G_x)(4):", nstr(8*(G_22(4)+G_x(4)), 20))
# U_3 = 8 G_2(s)+4 z''(s)
print("8 G_2(3):", nstr(8*G_2(3), 20))
# direct quadrant sums (partial, no tail)
U1q = U2q = U3q = mpf(0)
for k in range(1, 501):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_
        M = m0*log(m0)**2 + (k-m0)*log(k-m0)**2
        U1q += 4*power(k,-s-1)*L*log(k)
        U2q += 4*power(k,-s-2)*L*L
        U3q += 4*power(k,-s-1)*M
print("U1 partial direct:", nstr(U1q, 18))
print("U2 partial direct:", nstr(U2q, 18))
print("U3 partial direct:", nstr(U3q, 18))
'''
exec(code)
