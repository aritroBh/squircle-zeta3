code = r'''
from squircle import *
from c2 import G_2, G_l, G_x, G_22
from mpmath import mp, zeta, mpf, nstr, log, power, fsum, diff as mpdiff
setdps(60)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
s = 3
print("zeta''(2):", nstr(zetapp(2), 22))
print("zeta''(3):", nstr(zetapp(3), 22))
print("G_l(3):", nstr(G_l(3), 22))
print("G_22(4):", nstr(G_22(4), 22))
print("G_x(4):", nstr(G_x(4), 22))
print("G_2(3):", nstr(G_2(3), 22))
c2 = (2*s**2*zetapp(s-1) - 2*s**2*zetapp(s) - 8*s**2*G_l(s)
      + 4*s*(s+1)*(G_22(s+1)+G_x(s+1)) - 4*s*G_2(s))
print("c_2(3):", nstr(c2, 25))
print("direct route: 2.742211541922965596861203")
# check U_2 = 8(G_22+G_x)(4)+4 z''(3) vs direct
U2d = 8*(G_22(4)+G_x(4))+4*zetapp(3)
print("U_2 formula:", nstr(U2d, 20))
# direct U_2 partial
U2q = mpf(0)
for k in range(1, 801):
    for m0 in range(1, k):
        a_ = m0*log(m0); b_ = (k-m0)*log(k-m0)
        U2q += 4*power(k,-s-2)*(a_+b_)**2
print("U_2 partial direct:", nstr(U2q, 18))
'''
exec(code)
