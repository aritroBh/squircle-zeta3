# vjet.py -- VERIFIER's own code (independent of the solver's E_k formulas and of D15b's hand-expanded series).
# K_p = { F_p(x) <= 1 },  F_p(x) = sum_j c_j * exp(d*lam_j) * |<v_j, x>|^p,   d = p - 1.
# Boundary in polar form rho(th) = F_p(u_th)^(-1/p), perimeter = int sqrt(rho^2 + rho'^2) dth.
# All d-dependence is carried by generic truncated power series in d (automatic jets), so the d^k
# coefficient of the integrand is computed without any hand-derived expansion.
import mpmath as mp

NJ = 4  # jets d^0..d^3

def jmul(a, b):
    return [mp.fsum(a[i]*b[k-i] for i in range(k+1)) for k in range(NJ)]

def jexp(a):
    # exp of a jet: exp(a0) * exp(a - a0)
    x = [mp.mpf(0)] + list(a[1:])
    res = [mp.mpf(1)] + [mp.mpf(0)]*(NJ-1); pw = list(res)
    for j in range(1, NJ):
        pw = jmul(pw, x)
        res = [res[k] + pw[k]/mp.factorial(j) for k in range(NJ)]
    e0 = mp.exp(a[0])
    return [e0*r for r in res]

def jlog(a):
    x = [mp.mpf(0)] + [a[k]/a[0] for k in range(1, NJ)]
    res = [mp.mpf(0)]*NJ; pw = [mp.mpf(1)] + [mp.mpf(0)]*(NJ-1)
    for j in range(1, NJ):
        pw = jmul(pw, x)
        res = [res[k] + mp.mpf((-1)**(j+1))/j*pw[k] for k in range(NJ)]
    res[0] = mp.log(a[0])
    return res

def jinv(a):
    r = [1/a[0]] + [mp.mpf(0)]*(NJ-1)
    for k in range(1, NJ):
        r[k] = -mp.fsum(a[i]*r[k-i] for i in range(1, k+1))/a[0]
    return r

def jsqrt(a):
    l = jlog(a)
    return jexp([x/2 for x in l])

def jscal(c, a): return [c*x for x in a]
def jadd(a, b): return [x+y for x, y in zip(a, b)]

ONE_OVER_P = None
def one_over_p():
    return [mp.mpf((-1)**k) for k in range(NJ)]  # 1/(1+d)

def integrand_jet(th, vs, cs, lams):
    u = (mp.cos(th), mp.sin(th)); du = (-mp.sin(th), mp.cos(th))
    Fj = [mp.mpf(0)]*NJ; dFj = [mp.mpf(0)]*NJ
    for v, c, lam in zip(vs, cs, lams):
        g = v[0]*u[0] + v[1]*u[1]; gp = v[0]*du[0] + v[1]*du[1]
        s = abs(g); sg = 1 if g > 0 else -1
        # c * e^{d lam} * s^{1+d} = c*s*exp(d*(lam + log s))
        ls = mp.log(s) + lam
        term = jexp([mp.log(c*s), ls] + [mp.mpf(0)]*(NJ-2))
        Fj = jadd(Fj, term)
        # derivative: c e^{d lam} (1+d) s^d sgn gp
        t2 = jexp([mp.log(c), ls] + [mp.mpf(0)]*(NJ-2))
        t2 = jmul(t2, [mp.mpf(1), mp.mpf(1)] + [mp.mpf(0)]*(NJ-2))
        dFj = jadd(dFj, jscal(sg*gp, t2))
    ip = one_over_p()
    lF = jlog(Fj)
    rho = jexp(jscal(-1, jmul(ip, lF)))
    drho = jscal(-1, jmul(jmul(rho, ip), jmul(dFj, jinv(Fj))))
    return jsqrt(jadd(jmul(rho, rho), jmul(drho, drho)))

def kink_angles(vs):
    # directions th with <v,u_th> = 0, in [0, 2pi)
    out = []
    for v in vs:
        th0 = mp.atan2(v[1], v[0]) + mp.pi/2
        for k in range(2):
            out.append((th0 + k*mp.pi) % (2*mp.pi))
    return sorted(out)

def arc_coeffs(a, b, vs, cs, lams):
    """int_a^b of the d^k coefficient (k=0..3) of the polar integrand"""
    res = []
    for k in range(NJ):
        res.append(mp.quad(lambda th: integrand_jet(th, vs, cs, lams)[k], [a, b], method='tanh-sinh'))
    return res

def arc_coeffs_all(a, b, vs, cs, lams):
    # one quadrature per coefficient but re-using evaluations through a cache
    cache = {}
    def f(k):
        def g(th):
            key = th
            if key not in cache:
                cache[key] = integrand_jet(th, vs, cs, lams)
            return cache[key][k]
        return g
    return [mp.quad(f(k), [a, b], method='tanh-sinh') for k in range(NJ)]

def P_direct(p, vs, cs, lams):
    d = p - 1
    def f(th):
        u = (mp.cos(th), mp.sin(th)); du = (-mp.sin(th), mp.cos(th))
        F = mp.mpf(0); dF = mp.mpf(0)
        for v, c, lam in zip(vs, cs, lams):
            g = v[0]*u[0] + v[1]*u[1]; gp = v[0]*du[0] + v[1]*du[1]
            s = abs(g); sg = 1 if g > 0 else -1
            W = c*mp.exp(d*lam)
            F += W*s**p
            dF += W*p*s**(p-1)*sg*gp
        rho = F**(-1/p); drho = -(1/p)*F**(-1/p-1)*dF
        return mp.sqrt(rho**2 + drho**2)
    ks = kink_angles(vs)
    ks = [k for k in ks if k < mp.pi - mp.mpf(10)**(-mp.mp.dps//2)]
    pts = sorted(set([mp.mpf(0)] + ks + [mp.pi]))
    return 2*mp.quad(f, pts, method='tanh-sinh')
