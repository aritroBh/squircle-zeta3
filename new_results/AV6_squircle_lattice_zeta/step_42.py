code = r'''
"""Method C: boundary-layer/shell decomposition.
Z_p(s) = 4 zeta(s) (1 + 2^{-q}) + 8 T,   T = sum_{k>j>=1} (k^p+j^p)^{-q},  q = s/p.
T = sum_{k=2}^{K} k^{-s} inner(k)  +  tail_K,
inner(k) = sum_{j=1}^{k-1} (1+(j/k)^p)^{-q},
tail_K = I0 HT(K,s-1) + C0 HT(K,s)
         + sum_{r>=1} B_{2r}/(2r)! psi^{(2r-1)}(1) HT(K, s-1+2r)
         + sum_{j>=1} (-1)^j (q)_j/j! zeta(-jp) HT(K, s+jp),
I0 = 2F1(q, 1/p; 1+1/p; -1),  C0 = -(1+2^{-q})/2, psi(x) = (1+x^p)^{-q}.
"""
from mpmath import (mp, mpf, gamma, zeta, bernoulli, fac, pi, exp, log, sin,
                    power, fsum, hyper, loggamma)

def _psi_derivs_1(p, q, R):
    """psi^{(m)}(1), m=0..R, for psi(x)=(1+x^p)^{-q}, via series in t at x=1."""
    bk = [mpf(0)]*(R+1); b = mpf(1)
    for k in range(R+1):
        if k > 0: b *= (p-(k-1))/k
        bk[k] = b
    u = [bki/2 for bki in bk]; u[0] = mpf(0)
    W = [mpf(0)]*(R+1); W[0] = 1
    for n in range(R):
        acc = mpf(0)
        for j in range(1, n+2):
            acc += ((q-1)*j + n+1)*u[j]*W[n+1-j]
        W[n+1] = -acc/(n+1)
    twoq = power(2, -q)
    return [W[k]*twoq*fac(k) for k in range(R+1)]

def _poch(q, j):
    return gamma(q+j)/gamma(q)

def ZC(p, s, K=800, R=12, J=None, direct=True, info=False):
    p = mpf(p); s = mpf(s); q = s/p
    if J is None:
        J = max(2, int(2*R/float(p)) + 2)
    d1 = _psi_derivs_1(p, q, 2*R+1)
    I0 = hyper([q, 1/p], [1+1/p], -1)
    C0 = -(1+power(2, -q))/2
    def HT(K, sig):
        return zeta(sig) - fsum([power(m, -sig) for m in range(1, K+1)])
    P = mpf(0)
    if direct:
        for k in range(2, K+1):
            inner = mpf(0)
            kf = mpf(k)
            for j in range(1, k):
                inner += power(1+power(j/kf, p), -q)
            P += power(kf, -s)*inner
    tail = I0*HT(K, s-1) + C0*HT(K, s)
    for r in range(1, R+1):
        tail += bernoulli(2*r)/fac(2*r) * d1[2*r-1] * HT(K, s-1+2*r)
    cusp_terms = []
    for j in range(1, J+1):
        Cj = ((-1)**j) * _poch(q, j)/fac(j) * zeta(-j*p)
        cusp_terms.append(Cj*HT(K, s+j*p))
        tail += cusp_terms[-1]
    T = P + tail
    Z = 4*zeta(s)*(1+power(2, -q)) + 8*T
    if info:
        return Z, dict(P=P, tail=tail)
    return Z
'''
with open("methodC.py","w") as f: f.write(code)
test = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, dirichlet, mpf, nstr, fabs
setdps(60)
def beta(x): return dirichlet(x, [0,1,0,-1])
def rel(a,b):
    d = fabs(a-b); return d/fabs(b) if d else mpf('1e-999')
import time
t0=time.time()
# p=1 exact
print("p=1,s=4: ZC:", nstr(ZC(1,4,K=300), 40))
print("   exact 4z(3):", nstr(4*zeta(3), 40))
# p=2 exact
print("p=2,s=4: ZC:", nstr(ZC(2,4,K=300), 40))
print("   Epstein:", nstr(4*zeta(2)*beta(2), 40))
# p=50 vs ZA
t1=time.time()
print("p=50,s=4: ZC:", nstr(ZC(50,4,K=300,R=14), 40), f"{time.time()-t1:.0f}s")
print("        : ZA:", nstr(ZA(50,4), 40))
'''
exec(test)
with open("progress.log","a") as f: f.write("step_31: Method C implemented (methodC.py)\n")
