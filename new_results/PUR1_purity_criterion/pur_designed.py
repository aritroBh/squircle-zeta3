"""(f) designed square example and (g) designed rectangle example: identity, exp-small term,
expansion, convexity."""
import sympy as sp, mpmath as mp

mp.mp.dps = 50


def binom(t, n):
    return mp.binomial(t, n)


def Ttail(V, tau, nmax=400):
    """T(V) = int_V^inf e^v[(1+e^{-2v/tau})^tau - 1] dv = tau e^V sum_n binom(tau,n) e^{-2nV/tau}/(2n-tau)."""
    s = mp.mpf(0)
    for n in range(1, nmax):
        term = binom(tau, n) * mp.exp(-2 * n * V / tau) / (2 * n - tau)
        s += term
        if abs(term) < mp.mpf(10)**(-mp.mp.dps - 10):
            break
    return tau * mp.exp(V) * s


def tauG(tau):
    return 1 - mp.gamma(1 - tau / 2)**2 * mp.rgamma(1 - tau)


# ---------------- (f) square ----------------
Lam = lambda u: mp.log(2 + mp.sin(u) - mp.cos(u))


def h_sq(th, tau):
    j = mp.nint(th / (mp.pi / 2))
    u = abs(th - j * mp.pi / 2)
    return (mp.cos(u) + mp.sin(u)) * (1 + mp.exp(-2 * Lam(u) / tau))**tau


print("=== (f) designed square: int h = 16 + 2 tau Gamma(-tau/2)^2/Gamma(-tau) + E(tau) ===")
print("E(tau) = -16 tau sum_{n>=1} binom(tau,n) 4^{-n/tau}/(2n - tau)")
for tau in ['1', '0.7', '0.5', '0.3', '0.2', '0.1', '0.05', '0.02']:
    tau = mp.mpf(tau)
    pts = [mp.mpf(0)] + [mp.pi / 4 * mp.mpf(2)**(-k) for k in range(25, -1, -1)]
    I8 = 8 * mp.quad(lambda u: (mp.cos(u) + mp.sin(u)) * (1 + mp.exp(-2 * Lam(u) / tau))**tau, pts)
    E = -8 * Ttail(mp.log(2), tau)
    closed = 16 - 8 * mp.gamma(1 - tau / 2)**2 * mp.rgamma(1 - tau)
    print("tau=%-5s int h = %s ; closed+E - int h = %s ; E = %s" % (mp.nstr(tau, 3), mp.nstr(I8, 40), mp.nstr(closed + E - I8, 3), mp.nstr(E, 6)))
# one full-circle check, no symmetry reduction
tau = mp.mpf('0.3')
br = [k * mp.pi / 4 for k in range(9)]
Ifull = mp.quad(lambda th: h_sq(th, tau), br)
print("full-circle quadrature at tau=0.3, 8 pieces: %s  (diff to 8x reduced: %s)" % (
    mp.nstr(Ifull, 40), mp.nstr(Ifull - (16 - 8 * mp.gamma(1 - tau / 2)**2 * mp.rgamma(1 - tau) - 8 * Ttail(mp.log(2), tau)), 3)))

t = sp.Symbol('tau')
expo = sum(sp.zeta(k) * (sp.Rational(2)**(1 - k) - 1) * t**k / k for k in range(2, 9))
ser = sp.expand(sp.series(16 - 8 * sp.exp(expo), t, 0, 8).removeO())
zs = {k: sp.Symbol('zeta(%d)' % k) for k in range(2, 9)}
print("\nExpansion of 16 - 8 Gamma(1-tau/2)^2/Gamma(1-tau) through tau^7:")
for n in range(0, 8):
    cn = sp.simplify(ser.coeff(t, n))
    print("  tau^%d : %s   = %s" % (n, cn, mp.nstr(mp.mpf(sp.N(cn, 50)), 30)))
# numeric check of the series at small tau against the closed form
tau = mp.mpf('0.01')
sv = sum(mp.mpf(sp.N(ser.coeff(t, n), 60)) * tau**n for n in range(8))
print("check: closed(0.01) - series_7(0.01) = %s  (expect O(tau^8) ~ 1e-16)" % mp.nstr(16 - 8 * mp.gamma(1 - tau / 2)**2 * mp.rgamma(1 - tau) - sv, 5))

print("\nConvexity of the designed square body:")


def hpp_formula(u, tau):
    M, Mp = mp.cos(u) + mp.sin(u), mp.cos(u) - mp.sin(u)
    L = Lam(u)
    yy = 2 * L / tau
    sig = 1 / (1 + mp.exp(yy))
    f = (1 + mp.exp(-yy))**tau
    K = (1 - sig) / tau + sig
    return 2 * sig * f * M * mp.exp(-L) * (-3 * Mp + M**2 * mp.exp(-L) * (1 + 2 * K))


for tau in ['1', '0.8', '0.5', '0.2', '0.1', '0.03']:
    tau = mp.mpf(tau)
    hfun = lambda u: (mp.cos(u) + mp.sin(u)) * (1 + mp.exp(-2 * Lam(u) / tau))**tau
    # formula vs numerical second derivative at a few points
    errs = [abs(hpp_formula(u, tau) - (hfun(u) + mp.diff(hfun, u, 2))) for u in [mp.mpf('0.05'), mp.mpf('0.3'), mp.mpf('0.7')]]
    grid = [mp.pi / 4 * k / 400 for k in range(0, 401)]
    vals = [hpp_formula(u, tau) for u in grid]
    mn = min(vals)
    hp_left = mp.diff(hfun, mp.pi / 4)               # h'(pi/4^-) from the window around 0
    hp0 = mp.diff(hfun, mp.mpf(0), direction=1)      # one-sided derivative at 0
    print(" tau=%-4s max|formula - (h+h'') numeric| = %s ; min_{[0,pi/4]} (h+h'') = %s ; h'(0+) = %s ; jump of h' at pi/4 = %s" % (
        mp.nstr(tau, 3), mp.nstr(max(errs), 3), mp.nstr(mn, 6), mp.nstr(hp0, 3), mp.nstr(-2 * hp_left, 6)))

# ---------------- (g) designed rectangle ----------------
print("\n=== (g) designed rectangle [-A,A]x[-B,B], A=1, B=2 ===")
A, B = mp.mpf(1), mp.mpf(2)
L0 = lambda u: mp.log(1 + (B / A) * mp.sin(u) + (B / A)**2 * (1 - mp.cos(u)))
L1 = lambda w: mp.log(1 + (A / B) * mp.sin(w) + (A / B)**2 * (1 - mp.cos(w)))
us = mp.findroot(lambda u: L0(u) - L1(mp.pi / 2 - u), mp.mpf('0.5'))
V = L0(us)
c0, c1 = A**2 / B, B**2 / A
print("split angle u* = %s, V = %s ; c0 = A^2/B, c1 = B^2/A" % (mp.nstr(us, 30), mp.nstr(V, 30)))


def h_rect_q(th, tau):  # theta in [0, pi/2]
    M = A * mp.cos(th) + B * mp.sin(th)
    L = L0(th) if th <= us else L1(mp.pi / 2 - th)
    return M * (1 + mp.exp(-2 * L / tau))**tau


for tau in ['1', '0.5', '0.2', '0.1', '0.05']:
    tau = mp.mpf(tau)
    pts = sorted(set([mp.mpf(0)] + [us * mp.mpf(2)**(-k) for k in range(20, -1, -1)] + [mp.pi / 2 - (mp.pi / 2 - us) * mp.mpf(2)**(-k) for k in range(0, 21)] + [mp.pi / 2]))
    I = 4 * mp.quad(lambda th: h_rect_q(th, tau), pts)
    pred = 4 * (A + B) + 4 * (c0 + c1) * (tauG(tau) - Ttail(V, tau))
    print(" tau=%-4s int h = %s ; predicted - int h = %s" % (mp.nstr(tau, 3), mp.nstr(I, 40), mp.nstr(pred - I, 3)))
    # convexity grid
    hf = lambda th: h_rect_q(th, tau)
    g1 = [us * k / 120 for k in range(1, 120)]
    g2 = [us + (mp.pi / 2 - us) * k / 120 for k in range(1, 120)]
    mn = min(hf(x) + mp.diff(hf, x, 2) for x in g1 + g2)
    jump = mp.diff(hf, us, direction=1) - mp.diff(hf, us, direction=-1)
    d0 = mp.diff(hf, mp.mpf(0), direction=1)
    d90 = mp.diff(hf, mp.pi / 2, direction=-1)
    print("          min(h+h'') on grid = %s ; jump h'(u*+)-h'(u*-) = %s ; h'(0+) = %s ; h'(pi/2-) = %s" % (
        mp.nstr(mn, 6), mp.nstr(jump, 6), mp.nstr(d0, 3), mp.nstr(d90, 3)))
