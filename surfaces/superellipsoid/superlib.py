# superlib.py -- exact one-dimensional representation of the superellipsoid surface area.
#
#   S(p) = 8 A(eps),  eps = 1/p,
#   A(eps) = [3 G1^2 - 3 eps G1 H_2(eps) + eps^2 H_3(eps)] / Gamma(1+2 eps),   G1 = Gamma(1+eps),
#   H_m(eps) = (1/sqrt(pi)) int_0^oo h_eps(b)^m db,
#   h_eps(b) = int_0^oo (1 - exp(-x^(2-2eps)/b^2)) x^(eps-1) e^(-x) dx.
#
# h_eps is evaluated by its Mellin-Barnes residue series (convergent, small/moderate b)
#   h_eps(b) = Gamma(eps) - 1/(2-2eps) sum_{m>=0} (-1)^m Gamma(s_m) b^(2 s_m)/m!,  s_m=(m+eps)/(2-2eps)
# and its asymptotic series for large b
#   h_eps(b) ~ sum_{n>=1} (-1)^(n+1) Gamma(eps+(2-2eps)n) b^(-2n)/n!.
# At eps=0 exactly we use the closed form of k(b)=h_0(b) (see 01_kfun_Q0.py).
import mpmath as mp

def _bsw():
    return 3.1*mp.sqrt(mp.mp.dps)

def _asym_try(eps, b):
    """Asymptotic (n>=1) series in 1/b^2; returns (value, smallest |term|)."""
    b = mp.mpf(b); c = 2 - 2*eps
    s = 0; n = 1; prev = None; best = mp.inf
    tol = mp.mpf(10)**(-(mp.mp.dps + 8))
    lb = mp.log(b)
    while n < 200000:
        lt = mp.re(mp.loggamma(eps + c*n)) - mp.loggamma(n+1) - 2*n*lb
        if prev is not None and lt > prev and n > 3:
            break
        t = (-1)**(n+1)*mp.exp(mp.loggamma(eps + c*n) - mp.loggamma(n+1) - 2*n*lb)
        s += t
        at = abs(t)
        if at < best: best = at
        if at < tol*(1+abs(s)): break
        prev = lt; n += 1
    return s, best

def h_series(eps, b):
    b = mp.mpf(b)
    c = 2 - 2*eps
    lb2 = 2*mp.log(b)
    # scan term magnitudes to choose the working precision
    peak = mp.mpf(0); m = 0; ln10 = mp.log(10)
    while True:
        sm = (m + eps)/c
        lt = mp.re(mp.loggamma(sm) + sm*lb2) - mp.loggamma(m+1)
        if lt > peak: peak = lt
        if m > 10 and lt < peak - (mp.mp.dps + 30)*ln10 and lt < -(mp.mp.dps+30)*ln10:
            mmax = m; break
        m += 1
    extra = int(peak/ln10) + 20 + int(max(0, -mp.log10(abs(eps)+mp.mpf(10)**-300)))
    with mp.workdps(mp.mp.dps + extra):
        e = eps; c2 = 2 - 2*e; tot = 0
        for m in range(mmax+1):
            sm = (m + e)/c2
            tot += (-1)**m*mp.gamma(sm)*mp.exp(sm*2*mp.log(b))/mp.factorial(m)
        val = mp.gamma(e) - tot/c2
    return +val

def k0_closed(b):
    b = mp.mpf(b)
    with mp.workdps(mp.mp.dps + int(float(b*b)/4/2.3) + 10):
        v = -mp.log(b) - mp.euler/2 + mp.pi/2*mp.erfi(b/2) - b**2/4*mp.hyp2f2(1, 1, 1.5, 2, b**2/4)
    return +v

def h(eps, b):
    b = mp.mpf(b)
    if eps == 1:
        return -mp.expm1(-1/b**2)
    if eps == mp.mpf(1)/2:
        return mp.sqrt(mp.pi)*(1 - 1/mp.sqrt(1 + 1/b**2))
    if abs(mp.re(eps) - mp.mpf(1)/2) < mp.mpf(1)/8 and abs(b - 1) < mp.mpf(1)/2:
        e = eps
        return mp.quad(lambda x: -mp.expm1(-x**(2-2*e)/b**2)*x**(e-1)*mp.exp(-x), [0, b**(1/(1-e)), mp.inf])
    if mp.re(eps) >= mp.mpf(1)/2 or b > 1:
        val, best = _asym_try(eps, b)
        if best < mp.mpf(10)**(-(mp.mp.dps + 5)):
            return val
    if eps == 0:
        return k0_closed(b)
    return h_series(eps, b)

def H(m, eps):
    bs = _bsw()
    pts = [0, mp.mpf(1)/8, mp.mpf(1)/2, 1, 2, 4, 8, 16, bs, 2*bs, mp.inf] if bs > 16 else [0, mp.mpf(1)/8, mp.mpf(1)/2, 1, 2, 4, 8, 16, 32, mp.inf]
    return mp.quad(lambda b: h(eps, b)**m, pts)/mp.sqrt(mp.pi)

def A_oct(eps):
    eps = mp.mpf(eps)
    G1 = mp.gamma(1+eps)
    return (3*G1**2 - 3*eps*G1*H(2, eps) + eps**2*H(3, eps))/mp.gamma(1+2*eps)

def S_area(p):
    return 8*A_oct(1/mp.mpf(p))

def F2(eps):   # squircle: P(p)/4 = 2 - eps F(eps)
    return H(2, eps)/mp.gamma(1+eps)

def Qfun(eps):  # vertex function: Q(eps) = -H_3/Gamma(1+2eps)
    return -H(3, eps)/mp.gamma(1+2*eps)
