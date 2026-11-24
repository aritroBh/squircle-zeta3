code = r'''
"""AV6: squircle lattice zeta function  Z_p(s) = sum_{(m,n)!=(0,0)} (|m|^p+|n|^p)^{-s/p}.
Two independent evaluators:
  Method A: direct sum + one-dimensional Euler-Maclaurin strips (with algebraic/exponential
            Poisson-type tail corrections for sum_{n in Z} h(m,n) - integral).
  Method B: Mellin (theta) integral representation, split at t=delta; small-t part analytic
            via the asymptotic series of R_p(t); large-t part by termwise incomplete gamma
            or numerical quadrature.
Cross-checks: p=1 -> 4 zeta(s-1); p=2 -> 4 zeta(s/2) beta(s/2); p=inf -> 8 zeta(s-1).
"""
from mpmath import (mp, mpf, mpc, gamma, zeta, bernoulli, fac, pi, exp, log, sin, cos,
                    power, fsum, hyper, inf, quad, binomial, gammainc, loggamma, sqrt)

def setdps(d):
    mp.dps = d

def Aint(p, a):
    """int_0^inf (1+u^p)^{-a} du = (1/p) B(1/p, a-1/p)."""
    return gamma(1/p)*gamma(a-1/p)/(p*gamma(a))

def Jint(x, p, a):
    """int_x^inf (1+u^p)^{-a} du, stable in both regimes."""
    x = mpf(x)
    if x <= 1:
        return Aint(p,a) - x*hyper([a, 1/p],[1+1/p], -x**p)
    else:
        return x**(1-a*p)/(a*p-1)*hyper([a, a-1/p],[a-1/p+1], -x**(-p))

def taylor_h(m, y, p, a, N):
    """Series of g(t) = (m^p + (y+t)^p)^{-a} about t=0, to order N.
    Returns c[k] = g^{(k)}(y)/k!.  Works for real p>0, y>=1, m>=0."""
    y = mpf(y); m = mpf(m)
    v = [mpf(0)]*(N+2)
    yp = y**p
    for k in range(N+2):
        # binom(p,k) y^{p-k}: use loggamma form for general p, exact zero for integer p<k
        if float(p).is_integer() and k > int(p):
            bk = mpf(0)
        else:
            bk = exp(loggamma(p+1)-loggamma(k+1)-loggamma(p-k+1))
        v[k] = bk * yp / y**k
    c0 = v[0] + m**p
    v[0] = v[0] + m**p
    # normalize: v(t) = c0*(1+w), w = v/c0 - 1
    w = [vi/c0 for vi in v]
    w[0] -= 1
    # (1+w)^{-a}: recurrence (1+w) W' = -a w' W
    W = [mpf(0)]*(N+2)
    W[0] = 1
    for n in range(N+1):
        s_ = 0
        for j in range(1, n+2):
            s_ += ((a-1)*j + n+1) * w[j] * W[n+1-j]
        W[n+1] = -s_/(n+1)
    c0a = c0**(-a)
    return [W[k]*c0a for k in range(N+1)]

def harmtail(M, s):
    """sum_{m>M} m^{-s} = zeta(s) - H_M^{(s)}, computed stably."""
    return zeta(s) - fsum([power(m, -s) for m in range(1, M+1)])

def eps_terms(p, s, a, Imax):
    """Coefficients f_i of the asymptotic series of F(v)=int_R (1+|u|^p)^{-a} e^{2 pi i v u} du
       F(v) ~ sum_i f_i v^{-ip-1};  f_i = 2(-1)^{i+1} (a)_i sin(pi i p/2) Gamma(ip+1)/(i! (2pi)^{ip+1})."""
    fs = []
    for i in range(1, Imax+1):
        poch = gamma(a+i)/gamma(a)
        fi = 2*((-1)**(i+1))*poch*sin(pi*i*p/2)*gamma(i*p+1)/(fac(i)*(2*pi)**(i*p+1))
        fs.append(fi)
    return fs

def ZA(p, s, N=40, M=400, K=200, R=14, Icorr=12, info=False):
    """Method A: Euler-Maclaurin strips."""
    p = mpf(p); s = mpf(s)
    a = s/p
    Ap = Aint(p, a)
    def h(m, n):
        return power(power(m, p)+power(n, p), -a)
    # central quadrant sum
    Sc = fsum([h(m, n) for m in range(1, N+1) for n in range(1, N+1)])
    # strip E1' = sum_{m>N, n>=1} h(m,n)
    bern = [bernoulli(2*r)/fac(2*r) for r in range(R+1)]
    strip = mpf(0)
    for m in range(N+1, M+1):
        inner = fsum([h(m, n) for n in range(1, K+1)])
        # EM tail n>K
        ck = taylor_h(m, K, p, a, 2*R+1)
        tail = m*Jint(mpf(K)/m, p, a) - h(m, K)/2
        for r in range(1, R+1):
            tail -= bern[r]*ck[2*r-1]*fac(2*r-1)
        strip += inner + tail
    # tail m>M :  T1(m) = Ap*m^{1-s} - m^{-s}/2 + (1/2)eps(m), eps(m) ~ 2 sum_i f_i zeta(ip+1) m^{-s-ip}
    tailM = Ap*harmtail(M, s-1) - harmtail(M, s)/2
    fs = eps_terms(p, s, a, Icorr)
    for i, fi in enumerate(fs, start=1):
        tailM += fi*zeta(i*p+1)*harmtail(M, s+i*p)
    S = Sc + 2*strip
    Z = 4*zeta(s) + 4*S
    if info:
        return Z, dict(Sc=Sc, strip=strip, tailM=tailM)
    return Z

def r_series(p, s, Jmax):
    """r_j coefficients:  R_p(t) = f(t) - Gamma(1+1/p) t^{-1/p} + 1/2 ~ sum_{j>=1} r_j t^j,
       r_j = e_j zeta(jp+1), e_j = -2(-1)^j Gamma(jp+1) sin(pi j p/2)/(j! (2pi)^{jp+1})."""
    rs = []
    for j in range(1, Jmax+1):
        ej = -2*((-1)**j)*gamma(j*p+1)*sin(pi*j*p/2)/(fac(j)*(2*pi)**(j*p+1))
        rs.append(ej*zeta(j*p+1))
    return rs

def ZB(p, s, delta=mpf('0.02'), Jmax=80, D=100, info=False):
    """Method B: Mellin integral, small-delta split.
       Z = 1/Gamma(a) int_0^inf t^{a-1} (theta(t)^2-1) dt,  theta = 1+2 sum_{m>=1} e^{-t m^p}."""
    p = mpf(p); s = mpf(s); a = s/p
    ap = gamma(1+1/p)
    ga = gamma(a)
    # ---- small-t part on [0, delta]: analytic
    small = 4*ap*ap*delta**(a-2/p)/(a-2/p) - delta**a/a
    rs = r_series(p, s, Jmax)
    # truncate adaptively
    best = None; used = 0
    s1 = mpf(0)
    terms1 = []
    for j, rj in enumerate(rs, start=1):
        t1 = 8*ap*rj*delta**(a-1/p+j)/(a-1/p+j)
        terms1.append(abs(t1))
        s1 += t1
    # R^2 convolution
    s2 = mpf(0)
    conv = [mpf(0)]*(2*Jmax+1)
    for j in range(1, Jmax+1):
        for k in range(1, Jmax+1):
            conv[j+k] += rs[j-1]*rs[k-1]
    for j in range(2, 2*Jmax+1):
        s2 += 4*conv[j]*delta**(a+j)/(a+j)
    # find optimal truncation estimate of error (magnitude of smallest term in tail)
    small += s1 + s2
    # ---- large-t part on [delta, inf): termwise incomplete gamma
    # int_delta^inf t^{a-1}(theta^2-1)dt = 4 sum_{m,n>=1} b^{-a} Gamma(a, delta b) + 4 sum_{m>=1} m^{-s} Gamma(a, delta m^p)
    # where b = m^p+n^p.  cutoff b <= D/delta.
    L = int((D/delta)**(1/p)) + 2
    large1 = mpf(0)
    for m in range(1, L+1):
        mp_ = power(m, p)
        gm = gammainc(a, delta*mp_, inf)
        large1 += power(m, -s)*gm
    large1 *= 4
    large2 = mpf(0)
    mpows = [power(m, p) for m in range(L+1)]
    for m in range(1, L+1):
        mm = mpows[m]
        row = mpf(0)
        for n in range(1, L+1):
            b = mm + mpows[n]
            if b*delta > D*2:   # incomplete gamma negligible
                break
            row += power(b, -a)*gammainc(a, delta*b, inf)
        large2 += row
    large2 *= 4
    Z = (small + large1 + large2)/ga
    if info:
        return Z, dict(small=small, large1=large1, large2=large2, L=L)
    return Z
'''
with open("squircle.py","w") as f:
    f.write(code)
print("written", len(code))
with open("progress.log","a") as f:
    f.write("step_04: squircle.py module written (methods A and B)\n")
