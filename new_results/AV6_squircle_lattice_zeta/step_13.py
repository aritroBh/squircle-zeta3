code = r'''
"""AV6 squircle lattice zeta (v3, final core)."""
from mpmath import (mp, mpf, mpc, gamma, zeta, bernoulli, fac, pi, exp, log, sin, cos,
                    power, fsum, hyper, inf, quad, binomial, gammainc, loggamma, sqrt)

def setdps(d): mp.dps = d

def Aint(p, a):
    return gamma(1/p)*gamma(a-1/p)/(p*gamma(a))

def Jint(x, p, a):
    x = mpf(x)
    if x <= 1:
        return Aint(p,a) - x*hyper([a, 1/p],[1+1/p], -x**p)
    else:
        return x**(1-a*p)/(a*p-1)*hyper([a, a-1/p],[a-1/p+1], -x**(-p))

def taylor_h(m, y, p, a, N):
    y = mpf(y); m = mpf(m)
    v = [mpf(0)]*(N+2)
    yp = y**p
    pint = int(p) if float(p).is_integer() else -1
    for k in range(N+2):
        if k > pint >= 0:
            bk = mpf(0)
        else:
            bk = exp(loggamma(p+1)-loggamma(k+1)-loggamma(p-k+1))
        v[k] = bk * yp / y**k
    c0 = v[0] + m**p
    v[0] = c0
    w = [vi/c0 for vi in v]; w[0] -= 1
    W = [mpf(0)]*(N+2); W[0] = 1
    for n in range(N+1):
        acc = mpf(0)
        for j in range(1, n+2):
            acc += ((a-1)*j + n+1) * w[j] * W[n+1-j]
        W[n+1] = -acc/(n+1)
    c0a = c0**(-a)
    return [W[k]*c0a for k in range(N+1)]

def harmtail(M, s):
    return zeta(s) - fsum([power(m, -s) for m in range(1, M+1)])

def eps_terms(p, s, a, Imax):
    fs = []
    for i in range(1, Imax+1):
        poch = gamma(a+i)/gamma(a)
        fi = 2*((-1)**(i+1))*poch*sin(pi*i*p/2)*gamma(i*p+1)/(fac(i)*(2*pi)**(i*p+1))
        fs.append(fi)
    return fs

def ZA(p, s, N=40, M=None, K=250, R=16, Icorr=16, info=False):
    p = mpf(p); s = mpf(s); a = s/p
    D = mp.dps*mpf('2.35')
    if M is None:
        M = 300
        if abs(p-1) < mpf('0.49'):
            tau = abs(1-1/p)
            if tau > 0:
                M = max(M, int(1.3*float(D)/(2*pi*pi*tau))+20)
        if p > 2:
            M = max(M, int(1.3*float(D)/(2*pi*sin(pi/p)))+20)
    M = min(M, 200000)
    Ap = Aint(p, a)
    def h(m, n):
        return power(power(m, p)+power(n, p), -a)
    bern = [bernoulli(2*r)/fac(2*r) for r in range(R+1)]
    def T1(m):
        inner = fsum([h(m, n) for n in range(1, K+1)])
        ck = taylor_h(m, K, p, a, 2*R+1)
        tail = power(m, 1-s)*Jint(mpf(K)/m, p, a) - h(m, K)/2
        for r in range(1, R+1):
            tail -= bern[r]*ck[2*r-1]*fac(2*r-1)
        return inner + tail
    Sc = fsum([h(m, n) for m in range(1, N+1) for n in range(1, N+1)])
    strip = mpf(0); corner = mpf(0)
    for m in range(N+1, M+1):
        t1 = T1(m)
        strip += t1
        corner += t1 - fsum([h(m, n) for n in range(1, N+1)])
    tailM = Ap*harmtail(M, s-1) - harmtail(M, s)/2
    fs = eps_terms(p, s, a, Icorr)
    for i, fi in enumerate(fs, start=1):
        if fi != 0:
            tailM += fi*zeta(i*p+1)*harmtail(M, s+i*p)
    tc = tailM
    for n in range(1, N+1):
        cM = taylor_h(n, M, p, a, 2*R+1)
        em = power(n, 1-s)*Jint(mpf(M)/n, p, a) - h(M, n)/2
        for r in range(1, R+1):
            em -= bern[r]*cM[2*r-1]*fac(2*r-1)
        tc -= em
    strip += tailM
    corner += tc
    S = Sc + 2*strip - corner
    Z = 4*zeta(s) + 4*S
    if info: return Z, dict(Sc=Sc, strip=strip, corner=corner, M=M)
    return Z

def r_series(p, s, Jmax):
    rs = []
    for j in range(1, Jmax+1):
        ej = -2*((-1)**j)*gamma(j*p+1)*sin(pi*j*p/2)/(fac(j)*(2*pi)**(j*p+1))
        rs.append(ej*zeta(j*p+1))
    return rs

def ZB(p, s, delta=None, Jmax=200, info=False):
    """Method B: Mellin integral; small-t analytic, large-t quadrature."""
    p = mpf(p); s = mpf(s); a = s/p
    D = mp.dps*mpf('2.35')
    if delta is None:
        if p <= 1:
            delta = mpf('1.0')
        else:
            d1 = (2*pi/(1.3*D))**p
            kap = 2*pi*(1-1/p)*(2*pi/p)**(1/(p-1))*sin(pi/(2*(p-1)))
            d2 = (kap/(1.3*D))**(p-1)
            delta = min(mpf('0.02'), d1, d2)
    delta = mpf(delta)
    ap = gamma(1+1/p); ga = gamma(a)
    small = 4*ap*ap*delta**(a-2/p)/(a-2/p) - delta**a/a
    rs = r_series(p, s, Jmax)
    # linear R part with optimal truncation
    terms = [8*ap*rs[j-1]*delta**(a-1/p+j)/(a-1/p+j) for j in range(1, Jmax+1)]
    s1 = mpf(0); run = mpf(0); best = inf; jstar = 0
    for j, tj in enumerate(terms, start=1):
        run += tj
        if abs(tj) < best:
            best = abs(tj); jstar = j; s1 = run
    est1 = best
    # R^2 part truncated at same order
    conv = [mpf(0)]*(2*Jmax+1)
    for j in range(1, Jmax+1):
        rj = rs[j-1]
        if rj == 0: continue
        for k in range(1, Jmax+1):
            conv[j+k] += rj*rs[k-1]
    s2 = mpf(0)
    for j in range(2, min(2*jstar, 2*Jmax)+1):
        if conv[j] != 0:
            s2 += 4*conv[j]*delta**(a+j)/(a+j)
    small += s1 + s2
    # large-t quadrature
    def theta(t):
        f = mpf(0)
        m = 1
        while True:
            e = exp(-t*power(m, p))
            f += e
            if e < mpf(10)**(-mp.dps-10): break
            m += 1
        return 1+2*f
    def integrand(t):
        th = theta(t)
        return t**(a-1)*(th*th-1)
    large = quad(integrand, [delta, 1, inf])
    Z = (small + large)/ga
    if info: return Z, dict(small=small, large=large, delta=delta, est1=est1, jstar=jstar)
    return Z
'''
with open("squircle.py","w") as f: f.write(code)

test = r'''
from squircle import *
from mpmath import mp, zeta, dirichlet, mpf, nstr, fabs, sin, pi
import time
setdps(60)
def beta(x): return dirichlet(x, [0,1,0,-1])
def rel(a,b):
    d = fabs(a-b)
    return d/fabs(b) if d else mpf('1e-999')
t0=time.time()
Z_B, infB = ZB(2, 4, info=True)
print("p=2,s=4 ZB:", nstr(Z_B,50), " vs exact rel:", nstr(rel(Z_B, 4*zeta(2)*beta(2)),4))
print("   delta:", nstr(infB['delta'],4), " est1:", nstr(infB['est1'],4), " jstar:", infB['jstar'])
Z_B3, infB3 = ZB(3, 4, info=True)
Z_A3, infA3 = ZA(3, 4, info=True)
print("p=3,s=4: ZB:", nstr(Z_B3,50)); print("         ZA:", nstr(Z_A3,50))
print("   agree:", nstr(rel(Z_B3,Z_A3),4), " delta:", nstr(infB3['delta'],4), " est1:", nstr(infB3['est1'],4), f" M={infA3['M']}")
print("time", time.time()-t0)
'''
exec(test)
with open("progress.log","a") as f: f.write("step_13: v3 module, quad tail + adaptive delta/M; p=2,3 checks\n")
