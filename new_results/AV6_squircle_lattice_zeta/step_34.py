code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, power, fsum, hyper, exp, log, binomial, fac
import time
setdps(70)
p = mpf(50); s = mpf(4); eps = 1/p; q = s*eps  # q = s/p = a
a = q
t0=time.time()
# psi(x) = (1+x^p)^(-q); derivatives at x=1 via series: psi(1+t) = 2^-q (1+u)^(-q), u = ((1+t)^p-1)/2
def psi_derivs(R):
    # series of (1+t)^p
    bk = [mpf(0)]*(R+1)
    b = mpf(1)
    for k in range(R+1):
        if k>0: b *= (p-(k-1))/k
        bk[k] = b
    u = [bk[k]/2 for k in range(R+1)]; u[0] = 0
    # (1+u)^(-q) via recurrence (1+u)W' = -q u' W
    W = [mpf(0)]*(R+1); W[0]=1
    for n in range(R):
        acc = mpf(0)
        for j in range(1, n+2):
            acc += ((q-1)*j + n+1)*u[j]*W[n+1-j]
        W[n+1] = -acc/(n+1)
    twoq = power(2, -q)
    return [W[k]*twoq*fac(k) for k in range(R+1)]  # psi^{(k)}(1)
Rc = 13
dpsi = psi_derivs(Rc)
I0 = hyper([q, 1/p],[1+1/p], -1)  # int_0^1 psi
psi0 = mpf(1); psi1 = power(2,-q)
K = 4000
# partial sum
P = mpf(0)
for k in range(2, K+1):
    inner = mpf(0)
    kj = mpf(1)/k
    for j in range(1, k):
        x = j*kj
        inner += power(1+power(x,p), -q)
    P += power(k, -s)*inner
    if k % 1000 == 0: print("k", k, time.time()-t0)
# tail with EM in k
tail = I0*harmtail(K, s-1) - (psi0+psi1)/2*harmtail(K, s)
sign = 1
for r in range(1, Rc//2+1):
    # term: (-1)^{r+1} B_{2r}/(2r)! psi^{(2r-1)}(1) k^{-(2r-1)} ... standard EM: +B_{2r}/(2r)! (g^{(2r-1)}(k)) with sign pattern from EM
    pass
# EM: sum_{j<=k-1} g(j) = k int psi - (psi0+psi1)/2 + sum_{r>=1} B_{2r}/(2r)! (psi^{(2r-1)}(1) k^{1-2r} - psi^{(2r-1)}(0) k^{1-2r})
# psi^{(m)}(0)=0 for m < p (p=50 > 13) so:
bern = [mpf(0)] + [bernoulli(2*r)/fac(2*r) for r in range(1, Rc//2+1)]
for r in range(1, Rc//2+1):
    tail += bern[r]*dpsi[2*r-1]*harmtail(K, s-1+2*r)
T = P + tail
ZC = 4*zeta(s)*(1+power(2,-q)) + 8*T
print("Method C:", nstr(ZC, 45), f"({time.time()-t0:.0f}s)")
print("ZA      :", nstr(ZA(50,4), 45))
print("ZB      : 9.38291268702971115535336280772466783686741258")
'''
with open("step_26.py","w") as f: f.write(code)
exec(code)
with open("progress.log","a") as f: f.write("step_26: Method C boundary-layer arbitration p=50\n")
