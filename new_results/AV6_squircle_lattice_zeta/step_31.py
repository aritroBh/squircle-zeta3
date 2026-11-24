code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, power, fsum, bernoulli, fac
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
def h(m,n): return power(power(m,p)+power(n,p),-a)
m = 41; K = 250; R = 16
bern = [bernoulli(2*r)/fac(2*r) for r in range(R+1)]
ck = taylor_h(m, K, p, a, 2*R+1)
tail = power(m, 1-s)*Jint(mpf(K)/m, p, a) - h(m, K)/2
for r in range(1, R+1):
    tail -= bern[r]*ck[2*r-1]*fac(2*r-1)
seg = fsum([h(m,n) for n in range(K+1, 2001)])
tail2000 = power(m,1-s)*Jint(mpf(2000)/m, p, a)
print("tail_K    :", nstr(tail, 35))
print("seg+t2000 :", nstr(seg+tail2000, 35))
print("diff      :", nstr(tail-seg-tail2000, 6))
# check pieces of the EM tail
print("int part  :", nstr(power(m,1-s)*Jint(mpf(K)/m,p,a), 35))
print("-h/2      :", nstr(-h(m,K)/2, 12))
bterm = mpf(0)
for r in range(1, R+1):
    bterm -= bern[r]*ck[2*r-1]*fac(2*r-1)
print("bern terms:", nstr(bterm, 12))
# exact: int_K^inf g - sum_{n>K} g(n) should = g(K)/2 + sum B terms
lhs = power(m,1-s)*Jint(mpf(K)/m,p,a) - seg - tail2000
print("int - sum :", nstr(lhs, 35), " (should = g(K)/2 + Bterms =", nstr(h(m,K)/2 - bterm, 35), ")")
'''
exec(code)
