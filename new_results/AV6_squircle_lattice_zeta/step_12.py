code = r'''
"""AV6 squircle lattice zeta (v2)."""
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
    w = [vi/c0 for vi in v]
    w[0] -= 1
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

def ZA(p, s, N=40, M=400, K=250, R=16, Icorr=16, info=False):
    """Method A: central square + 2 strips - corner, Euler-Maclaurin tails."""
    p = mpf(p); s = mpf(s)
    a = s/p
    Ap = Aint(p, a)
    def h(m, n):
        return power(power(m, p)+power(n, p), -a)
    bern = [bernoulli(2*r)/fac(2*r) for r in range(R+1)]
    def T1(m):
        """sum_{n>=1} h(m,n) = inner + EM tail n>K"""
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
    # tail m>M of sum T1(m)
    tailM = Ap*harmtail(M, s-1) - harmtail(M, s)/2
    fs = eps_terms(p, s, a, Icorr)
    for i, fi in enumerate(fs, start=1):
        tailM += fi*zeta(i*p+1)*harmtail(M, s+i*p)
    # tailcorner = sum_{m>M} [T1(m) - sum_{n<=N} h(m,n)]
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
    if info:
        return Z, dict(Sc=Sc, strip=strip, corner=corner, tailM=tailM)
    return Z

def r_series(p, s, Jmax):
    rs = []
    for j in range(1, Jmax+1):
        ej = -2*((-1)**j)*gamma(j*p+1)*sin(pi*j*p/2)/(fac(j)*(2*pi)**(j*p+1))
        rs.append(ej*zeta(j*p+1))
    return rs

def ZB(p, s, delta=mpf('0.02'), Jmax=80, D=100, info=False):
    """Method B: Mellin integral split at delta."""
    p = mpf(p); s = mpf(s); a = s/p
    ap = gamma(1+1/p)
    ga = gamma(a)
    small = 4*ap*ap*delta**(a-2/p)/(a-2/p) - delta**a/a
    rs = r_series(p, s, Jmax)
    s1 = mpf(0)
    for j, rj in enumerate(rs, start=1):
        s1 += 8*ap*rj*delta**(a-1/p+j)/(a-1/p+j)
    s2 = mpf(0)
    conv = [mpf(0)]*(2*Jmax+1)
    for j in range(1, Jmax+1):
        rj = rs[j-1]
        if rj == 0: continue
        for k in range(1, Jmax+1):
            conv[j+k] += rj*rs[k-1]
    for j in range(2, 2*Jmax+1):
        if conv[j] != 0:
            s2 += 4*conv[j]*delta**(a+j)/(a+j)
    small += s1 + s2
    L = int((D/delta)**(1/p)) + 2
    large1 = mpf(0)
    mpows = [power(m, p) for m in range(L+1)]
    for m in range(1, L+1):
        large1 += power(m, -s)*gammainc(a, delta*mpows[m], inf)
    large1 *= 4
    large2 = mpf(0)
    for m in range(1, L+1):
        row = mpf(0)
        mm = mpows[m]
        for n in range(1, L+1):
            b = mm + mpows[n]
            if b*delta > 2*D: break
            row += power(b, -a)*gammainc(a, delta*b, inf)
        large2 += row
    large2 *= 4
    Z = (small + large1 + large2)/ga
    if info:
        return Z, dict(small=small, large1=large1, large2=large2, L=L)
    return Z
'''
with open("squircle.py","w") as f: f.write(code)
test = r'''
from squircle import *
from mpmath import mp, zeta, dirichlet, mpf, nstr, fabs
setdps(60)
def beta(x): return dirichlet(x, [0,1,0,-1])
exact = 4*zeta(2)*beta(2)
Z_A = ZA(2, 4, N=30, M=300, K=150, R=12)
print("ZA   :", nstr(Z_A, 50))
print("exact:", nstr(exact, 50))
print("rel err:", nstr(fabs(Z_A-exact)/exact, 5))
'''
exec(test)
with open("progress.log","a") as f: f.write("step_12: Method A v2 with corner, retest p=2\n")
