code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, power, fsum, bernoulli, fac
setdps(70)
p = mpf(50); s = mpf(4); a = s/p
def h(m,n): return power(power(m,p)+power(n,p),-a)
def T1(m, K, R):
    bern = [bernoulli(2*r)/fac(2*r) for r in range(R+1)]
    inner = fsum([h(m, n) for n in range(1, K+1)])
    ck = taylor_h(m, K, p, a, 2*R+1)
    tail = power(m, 1-s)*Jint(mpf(K)/m, p, a) - h(m, K)/2
    for r in range(1, R+1):
        tail -= bern[r]*ck[2*r-1]*fac(2*r-1)
    return inner + tail, tail
m = 41
vb = fsum([h(m,n) for n in range(1,2001)]) + power(m,1-s)*Jint(mpf(2000)/m,p,a)
for K,R in [(250,16),(300,20),(400,20),(250,24)]:
    v1, tail = T1(m,K,R)
    print(f"K={K} R={R}: T1={nstr(v1,30)} diff={nstr(v1-vb,4)} tail={nstr(tail,12)}")
# check the derivative values c_k vs numerical diff at K=250
K=250
ck = taylor_h(m, K, p, a, 6)
g = lambda t: (m**p + (K+t)**p)**(-a)
from mpmath import diff as mpdiff
for k in [1,3,5]:
    print(k, nstr(ck[k]*fac(k), 15), nstr(mpdiff(g, 0, n=k), 15))
'''
exec(code)
