code = r'''
from squircle import *
from c2 import G_l
from mpmath import mp, zeta, mpf, nstr, log, power, fsum, diff as mpdiff
setdps(60)
def zetap(a): return mpdiff(lambda x: zeta(x), a)
s = 3
K = 5000
P = mpf(0); acc = mpf(0)
for k in range(1, K+1):
    if k >= 2: P += power(k, -s-1)*log(k)*acc
    acc += k*log(k)
# tail: sum_{m<=K} m ln m HTln(K,s+1) + Hurwitz part (same machinery as c2)
def HTln(sig):
    return -zetap(sig) - fsum([power(m,-sig)*log(m) for m in range(1, K+1)])
def HTln2(sig):
    return zetapp(sig) - fsum([power(m,-sig)*log(m)**2 for m in range(1, K+1)])
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
tail = HTln(s+1)*acc + (1/s)*HTln2(s-1) + (1/s**2)*HTln(s-1) + mpf(1)/2*HTln2(s)
print("P:", nstr(P, 20))
print("P+tail:", nstr(P+tail, 20))
print("G_l(3) from code:", nstr(G_l(3), 20))
'''
exec(code)
