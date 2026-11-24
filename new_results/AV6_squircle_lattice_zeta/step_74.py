code = r'''
from squircle import setdps
from Uq_direct import _lnq_derivs_1, stirling1u
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, exp, fsum, pi, diff, bernoulli, fac
setdps(60)
s = 4
p = 128; K = 768
def inner(k, p):
    kf = mpf(k)
    return fsum([log(1+power(1-mpf(d)/kf, p)) for d in range(1, k)])
# Uq_direct internals
I0q = fsum([(-1)**(n+1)/(n*(n*p+1)) for n in range(1, 400)])
d1 = _lnq_derivs_1(p, 1, 21)
# check derivatives against numerical
g = lambda x: log(1+x**p)
for m in [1, 3, 5]:
    print(f"psi^({m})(1): series={nstr(d1[m], 15)}  numeric={nstr(diff(g, 1, n=m), 15)}")
print("I0q:", nstr(I0q, 20))
# P partial sum both ways
P = mpf(0)
for k in range(2, K+1):
    P += power(k, -s)*inner(k, p)
print("P(K=768):", nstr(P, 25))
# tail
def HT(K, sig):
    return zeta(sig) - fsum([power(m, -sig) for m in range(1, K+1)])
tail = I0q*HT(K, s-1) - log(2)/2*HT(K, s)
for r in range(1, 11):
    tail += bernoulli(2*r)/fac(2*r) * d1[2*r-1] * HT(K, s-1+2*r)
print("tail:", nstr(tail, 25))
print("P+tail:", nstr(P+tail, 25))
# inner(k) EM check for k=768: EM expansion vs direct
em = K*I0q - log(2)/2 + sum(bernoulli(2*r)/fac(2*r)*d1[2*r-1]*K**(1-2*r) for r in range(1,11))
print("inner(768) direct:", nstr(inner(768,p), 20))
print("inner(768) EM    :", nstr(em, 20))
'''
exec(code)
