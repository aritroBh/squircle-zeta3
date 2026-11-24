code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum, diff as mpdiff
import mpmath
setdps(60)
s = 3
K = 1500
# U_1 = sum (|m|+|n|)^{-s-1} L ln S.  Quadrant part: 8 sum k^{-s-1} sigma(k) ln k
U1q = mpf(0)
acc = mpf(0)
for k in range(1, K+1):
    if k >= 2:
        U1q += 8*power(k, -s-1)*log(k)*acc
    acc += k*log(k)
# axes: 4 zeta''(s) = 4 sum k^{-s} ln^2 k
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
axes = 4*(zetapp(s) - fsum([power(m,-s)*log(m)**2 for m in range(1, K+1)]))
# tail: sigma(k) ln k ~ (k^2/2 ln^2 k - k^2/2 ln k + ...) * k^{-s-1}: use EM-tail:
# tail = sum_{m<=K} m ln m * HTln(K,s+1) + Hurwitz part
def HTln(sig):
    return -mpdiff(lambda x: zeta(x), sig) - fsum([power(m,-sig)*log(m) for m in range(1, K+1)])
def HTln2(sig):
    return zetapp(sig) - fsum([power(m,-sig)*log(m)**2 for m in range(1, K+1)])
tail = 8*( HTln(s+1)*acc + (1/s)*HTln2(s-1) + (1/s**2)*HTln(s-1) + mpf(1)/2*HTln2(s) )
# the 8x is because G_l tail was computed for sum k^{-s-1} ln k sigma(k) and U_1 = 8 G_l
print("U1 quadrant:", nstr(U1q, 20))
print("U1 axes    :", nstr(axes, 15))
print("U1 tail    :", nstr(tail, 15))
U1 = U1q + axes + tail
print("U_1 direct:", nstr(U1, 25))
print("needed for identity: 6.25064 ; G-formula: 6.26274")
'''
exec(code)
