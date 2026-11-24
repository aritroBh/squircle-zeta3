
"""AV6 squircle lattice zeta (v3, final core)."""
from mpmath import (mp, mpf, mpc, gamma, zeta, bernoulli, fac, pi, exp, log, sin, cos,
                    power, fsum, hyper, inf, quad, binomial, gammainc, loggamma, sqrt)

def setdps(d): mp.dps = d

def Aint(p, a):
    return gamma(1/p)*gamma(a-1/p)/(p*gamma(a))

def Jint(x, p, a):
    x = mpf(x)
    if x <= 1:
        return (Aint(p,a) - x*hyper([a, 1/p],[1+1/p], -x**p)).real
    else:
        return (x**(1-a*p)/(a*p-1)*hyper([a, a-1/p],[a-1/p+1], -x**(-p))).real

def taylor_h(m, y, p, a, N):
    y = mpf(y); m = mpf(m)
    v = [mpf(0)]*(N+2)
    yp = y**p
    pint = int(p) if float(p).is_integer() else -1
    bk = mpf(1)
    for k in range(N+2):
        if k > pint >= 0:
            bk = mpf(0)
        else:
            if k > 0:
                bk *= (p-(k-1))/k
            v[k] = bk * yp / y**k
        if k > pint >= 0:
            v[k] = mpf(0)
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

def ZB(p, s, delta=None, Jmax=None, info=False):
    if Jmax is None:
        Jmax = 160 if mpf(p) <= 1 else 200
    """Method B: Mellin integral; small-t analytic, large-t quadrature."""
    p = mpf(p); s = mpf(s); a = s/p
    D = mp.dps*mpf('2.35')
    if delta is None:
        if p <= 1:
            delta = mpf('4.0')
        elif p < 1.5:
            delta = min(mpf('0.02'), (2*pi/(1.3*D))**p)
        else:
            d1 = (2*pi/(1.3*D))**p
            kap = 2*pi*(1-1/p)*(2*pi/p)**(1/(p-1))*sin(pi/(2*(p-1)))
            d2 = (abs(kap)/(1.3*D))**(p-1)
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
    mpows_cache = {}
    def theta(t):
        f = mpf(0)
        m = 1
        cutoff = mpf(10)**(-mp.dps-10)
        while True:
            mpw = mpows_cache.get(m)
            if mpw is None:
                mpw = power(m, p); mpows_cache[m] = mpw
            e = exp(-t*mpw)
            f += e
            if e < cutoff: break
            m += 1
        return 1+2*f
    def F(t):
        th = theta(t)
        return th*th-1
    if delta < 1:
        ld = float(log(delta))
        if ld < -30:
            # chunked tanh-sinh on the log-interval (accuracy on long intervals)
            pts = [ld]
            while pts[-1] < -25:
                pts.append(pts[-1]+25)
            pts.append(0.0)
            I1 = mpf(0)
            for i in range(len(pts)-1):
                I1 += quad(lambda x: exp(a*x)*F(exp(x)), [pts[i], pts[i+1]])
        else:
            I1 = quad(lambda x: exp(a*x)*F(exp(x)), [ld, 0.0])
        I2 = quad(lambda t: t**(a-1)*F(t), [1, inf])
        large = I1 + I2
    else:
        large = quad(lambda t: t**(a-1)*F(t), [delta, inf])
    Z = (small + large)/ga
    if info: return Z, dict(small=small, large=large, delta=delta, est1=est1, jstar=jstar)
    return Z
