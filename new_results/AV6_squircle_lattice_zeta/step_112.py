code = r'''
"""G_2^{(2)}(s) = sum_k k^{-s-1} sum_{m<k} m^2 ln^2 m."""
from mpmath import mp, mpf, zeta, gamma, fac, binomial, log, power, nstr, fsum, bernoulli, diff as mpdiff

def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)

def HT(K, sig):
    return zeta(sig) - fsum([power(m,-sig) for m in range(1, K+1)])
def HTln2(K, sig):
    return zetapp(sig) - fsum([power(m,-sig)*log(m)**2 for m in range(1, K+1)])

def G_22(s, K=2000, R=10):
    P = mpf(0)
    acc = mpf(0)
    for k in range(1, K+1):
        if k >= 2: P += power(k, -s-1)*acc
        acc += k**2*log(k)**2
    tail = HT(K, s+1)*acc
    tail += (1/s)*HTln2(K, s-2) - mpf(1)/2*HTln2(K, s-1)
    for r in range(1, R+1):
        poch = gamma(s+1+2*r-1)/gamma(s+1)
        tail -= bernoulli(2*r)/fac(2*r)*poch*HTln2(K, s+2*r-2)
    return P + tail
'''
src = open("c2.py").read()
src += code
open("c2.py","w").write(src)
test = r'''
from c2 import *
from squircle import setdps
from mpmath import mp, nstr, zeta, diff as mpdiff
setdps(60)
def zetapp(a): return mpdiff(lambda x: zeta(x), a, 2)
def c2_exact(s):
    return (2*s**2*zetapp(s-1) - 2*s**2*zetapp(s) - 8*s**2*G_l(s)
            + 4*s*(s+1)*(G_22(s+1)+G_x(s+1)) - 4*s*G_2(s))
for s in [3,4,5,6]:
    print(f"s={s}: c_2={nstr(c2_exact(s), 30)}")
print("fits: 2.73718981077298834655 ; 1.573849628813923486848 ; 1.003456318418583651585 ; 0.646106636259136217485")
'''
exec(test)
