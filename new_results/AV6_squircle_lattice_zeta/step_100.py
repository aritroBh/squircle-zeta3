code = r'''
from c2 import *
from squircle import setdps
from mpmath import mp, nstr, zeta, log, power, fsum, mpf, diff as mpdiff
setdps(60)
s = 4
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
# direct shell sums for U_1 = sum T w lnS, U_2 = sum T w^2, U_3 = sum T (M/S)
K = 1200
U1 = U2 = U3 = mpf(0)
for k in range(1, K+1):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        L = a_+b_; w = L/k
        M = m0*log(m0)**2 + (k-m0)*log(k-m0)**2
        U1 += 2*power(k,-s-1)*L*log(k)
        U2 += 2*power(k,-s-2)*L*L
        U3 += 2*power(k,-s-1)*M
# axes contributions: |m|=k,n=0: L = k ln k, w = ln k, M = k ln^2 k
# U1: sum (|m|+|n|)^{-s-1} L ln S : axes: k^{-s-1} * k ln k * ln k * 2 = 2 k^{-s} ln^2 k
U1 += 2*(zetapp(s) - fsum([power(m,-s)*log(m)**2 for m in range(1,K+1)]))
U2 += 2*(zetapp(s) - fsum([power(m,-s)*log(m)**2 for m in range(1,K+1)]))
U3 += 2*(zetapp(s) - fsum([power(m,-s)*log(m)**2 for m in range(1,K+1)]))
print("U1 direct:", nstr(U1, 20), " (tails omitted)")
print("U2 direct:", nstr(U2, 20))
print("U3 direct:", nstr(U3, 20))
# G-based: U_1 = 8 G_l + 4 zeta''(s)
print("8 G_l + 4 z''(s):", nstr(8*G_l(s)+4*zetapp(s), 20))
print("8 G_2(s+1)+8 G_x(s+1)+2 z''(s):", nstr(8*G_2(s+1)+8*G_x(s+1)+2*zetapp(s), 20))
print("8 G_2(s) + 4 z''(s):", nstr(8*G_2(s)+4*zetapp(s), 20))
'''
exec(code)
