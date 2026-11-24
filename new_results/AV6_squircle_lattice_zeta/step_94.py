code = r'''
"""G(s) = sum_{k>=2} k^{-s-1} sigma(k), sigma(k) = sum_{m<k} m ln m.  Exact direct + EM tail."""
from mpmath import mp, mpf, zeta, gamma, fac, binomial, log, power, nstr, fsum, bernoulli, glaiser, diff as mpdiff

def zetap(a): return mpdiff(lambda x: zeta(x), a)

def G_sigma(s, K=2000):
    lnA = log(glaiser)
    P = mpf(0)
    for k in range(2, K+1):
        P += power(k, -s-1)*fsum([m0*log(m0) for m0 in range(1, k)])
    def HT(sig):
        return zeta(sig) - fsum([power(m,-sig) for m in range(1, K+1)])
    def HTln(sig):
        return -zetap(sig) - fsum([power(m,-sig)*log(m) for m in range(1, K+1)])
    T = mpf(0)
    T += mpf(1)/2*HTln(s-1) - mpf(1)/2*HTln(s) + mpf(1)/12*HTln(s+1)
    T += -mpf(1)/2*HT(s-1) + mpf(1)/2*HT(s) - mpf(1)/4*HT(s+1)
    T += lnA*HT(s+1)
    for i in range(1, 41):
        T += -(mpf(1)/(2*i))*HT(s-1+i) + (mpf(1)/(2*i))*HT(s+i) - (mpf(1)/(12*i))*HT(s+1+i)
    # R(k): lnH(k-1) correction: sum_j g_j (k-1)^{-2j}, g_j = B_{2j+2}/(2j(2j+2))
    for j in range(1, 13):
        cj = bernoulli(2*j+2)/(2*j*(2*j+2))
        for t in range(0, 41):
            T += cj*binomial(2*j+t-1, t)*HT(s+1+2*j+t)
    return P + T

def c1_exact(s, K=2000):
    return 4*s*(zetap(s) - zetap(s-1)) - 8*s*G_sigma(s, K)
'''
with open("c1.py","w") as f: f.write(code)
test = r'''
from c1 import *
from squircle import setdps
from mpmath import mp, nstr, zeta, diff as mpdiff
setdps(70)
def zetap(a): return mpdiff(lambda x: zeta(x), a)
for s in [3,4,5,6]:
    g = G_sigma(s, K=1200)
    c1 = 4*s*(zetap(s)-zetap(s-1)) - 8*s*g
    print(f"s={s}: G(s)={nstr(g, 25)}")
    print(f"      c_1 = 4s(z'-z') - 8s G = {nstr(c1, 40)}")
'''
exec(test)
