code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs, log, power, fsum, binomial, bernoulli, fac, pi
setdps(50)
s = mpf(4)
def U1_direct(p, K=3000):
    # sum_{k>j>=1} k^{-s} ln(1+(j/k)^p), inner by log; tail estimated by integral approx
    tot = mpf(0)
    for k in range(2, K+1):
        inner = mpf(0)
        kf = mpf(k)
        for j in range(1, k):
            inner += log(1+power(j/kf, p))
        tot += power(kf, -s)*inner
    # tail: inner(k) ~ k*I0q - ln2/2 ; I0q = int_0^1 ln(1+x^p)dx = sum (-1)^{n+1}/(n(n p+1))
    I0q = fsum([(-1)**(n+1)/(n*(n*p+1)) for n in range(1, 400)])
    tail = I0q*harmtail(K, s-1) - log(2)/2*harmtail(K, s)
    return tot + tail
for p_ in [5, 10, 20]:
    u = U1_direct(p_)
    predA = zeta(s-1)*pi**2/(12*p_)   # eps*pi^2/12*zeta(s-1)
    print(f"p={p_}: U1={nstr(u,20)}   predA={nstr(predA,12)}  ratio={nstr(u/predA,8)}")
'''
exec(code)
