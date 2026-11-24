code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, power, fsum, bernoulli, fac, quad, inf
setdps(60)
# exact EM check for g(y) = y^-4 at K=250, R=16
K = 250; R = 16
bern = [bernoulli(2*r)/fac(2*r) for r in range(R+1)]
g = lambda y: mpf(y)**(-4)
exact = zeta(4) - fsum([mpf(n)**(-4) for n in range(1, K+1)])  # sum_{n>K} n^-4
tail = mpf(1)/(3*K**3) - g(K)/2
for r in range(1, R+1):
    dk = ((-1)**(2*r-1)) * fac(2*r+2)/fac(3) * mpf(K)**(-4-(2*r-1))  # g^{(2r-1)}(K) = (-1)^k (4)_k K^{-4-k}
    tail -= bern[r]*dk
print("exact sum_{n>K} n^-4:", nstr(exact, 30))
print("EM tail               :", nstr(tail, 30))
print("diff                  :", nstr(tail-exact, 6))
# so the convention is right for pure power. Now the real g at p=50, m=41:
p = mpf(50); a = mpf(4)/p
def h(n): return power(power(41,p)+power(n,p),-a)
exact2 = fsum([h(n) for n in range(K+1, 4001)]) + 41**(-3)*Jint(mpf(4000)/41, p, a)
ck = taylor_h(41, K, p, a, 2*R+1)
tail2 = 41**(-3)*Jint(mpf(K)/41, p, a) - h(K)/2
for r in range(1, R+1):
    tail2 -= bern[r]*ck[2*r-1]*fac(2*r-1)
print("exact2:", nstr(exact2, 30))
print("tail2 :", nstr(tail2, 30))
print("diff  :", nstr(tail2-exact2, 6))
# what does the true derivative series look like? compare c1 exact vs taylor_h
print("c1:", nstr(ck[1], 20))
c1e = -a*p*mpf(K)**(p-1)*(41**p+mpf(K)**p)**(-a-1)
print("c1e:", nstr(c1e, 20))
'''
exec(code)
