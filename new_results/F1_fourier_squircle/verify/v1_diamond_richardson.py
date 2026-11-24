"""
VERIFIER V1 (independent of the solver's s-representation and Watson-lemma coefficients).
Diamond end: a_k(tau) - A_k computed by DIRECT theta-quadrature of
    (8/pi) int_0^{pi/4} cos(th) cos(4k th) [ (1+tan(th)^q)^tau - 1 ] dth ,   q = 1/tau,
(using expm1/log1p, so no cancellation against A_k), at tau_j = tau0 * 2^-j, then Neville/Richardson
extrapolation to tau -> 0 of
    f(tau) = (-1)^k pi (a_k - A_k) / (sqrt2 tau^2)                 -> zeta(2)            (claim)
    g(tau) = (f - zeta2)/tau                                       -> zeta(3)            (claim)
    h(tau) = (g - zeta3)/tau                                       -> -(33/16+7k^2) zeta(4)
    i(tau) = (h + (33/16+7k^2) zeta4)/tau                          -> -[(3+16k^2) z5 + (7+16k^2) z2 z3]/4
Each step uses the previous (already independently confirmed) constant.
"""
import mpmath as mp, sys, json

mp.mp.dps = 90
pi = mp.pi; s2 = mp.sqrt(2)

def a_minus_A(k, tau):
    q = 1 / tau
    def f(th):
        t = mp.tan(th)
        x = mp.exp(q * mp.log(t))
        return mp.cos(th) * mp.cos(4 * k * th) * mp.expm1(tau * mp.log1p(x))
    top = pi / 4
    ds = [0, mp.mpf(1)/8, mp.mpf(1)/4, mp.mpf(1)/2, 1, 2, 4, 8, 16, 32, 64, 128]
    pts = [top - d * tau for d in ds if d * tau < top]
    pts = sorted(set([mp.mpf(0)] + pts))
    return 8 / pi * mp.quad(f, pts)

def neville0(xs, ys):
    """value at 0 of the interpolating polynomial through (xs, ys)"""
    n = len(xs); P = list(ys)
    for m in range(1, n):
        for i in range(n - m):
            P[i] = (xs[i + m] * P[i] - xs[i] * P[i + 1]) / (xs[i + m] - xs[i])
    return P[0]

z2, z3, z4, z5 = mp.zeta(2), mp.zeta(3), mp.zeta(4), mp.zeta(5)
ks = [int(a) for a in sys.argv[1:]] or [0, 1, 2, 3]
J = 16
results = {}
for k in ks:
    tau0 = mp.mpf('0.02') / (1 + k)
    taus = [tau0 / mp.mpf(2) ** j for j in range(J)]
    vals = [(-1) ** k * pi * a_minus_A(k, t) / s2 for t in taus]
    f = [v / t ** 2 for v, t in zip(vals, taus)]
    c2 = neville0(taus, f); c2b = neville0(taus[1:], f[1:])
    g = [(fv - z2) / t for fv, t in zip(f, taus)]
    c3 = neville0(taus, g); c3b = neville0(taus[1:], g[1:])
    c4claim = -(mp.mpf(33) / 16 + 7 * k * k) * z4
    hh = [(gv - z3) / t for gv, t in zip(g, taus)]
    c4 = neville0(taus, hh); c4b = neville0(taus[1:], hh[1:])
    c5claim = -((3 + 16 * k * k) * z5 + (7 + 16 * k * k) * z2 * z3) / 4
    ii = [(hv - c4claim) / t for hv, t in zip(hh, taus)]
    c5 = neville0(taus, ii); c5b = neville0(taus[1:], ii[1:])
    print(f"k={k} (tau0={mp.nstr(tau0,3)}, {J} levels down to {mp.nstr(taus[-1],3)})")
    print(f"   c2' - zeta(2)                 = {mp.nstr(c2 - z2, 3)}   (extrapolation stability {mp.nstr(c2 - c2b, 3)})")
    print(f"   c3' - zeta(3)                 = {mp.nstr(c3 - z3, 3)}   (stability {mp.nstr(c3 - c3b, 3)})")
    print(f"   c4' - (-(33/16+7k^2) zeta4)   = {mp.nstr(c4 - c4claim, 3)}   (stability {mp.nstr(c4 - c4b, 3)})   c4'={mp.nstr(c4,25)}")
    print(f"   c5' - claim                   = {mp.nstr(c5 - c5claim, 3)}   (stability {mp.nstr(c5 - c5b, 3)})   c5'={mp.nstr(c5,25)}")
    sys.stdout.flush()
    results[k] = [mp.nstr(c2 - z2, 3), mp.nstr(c3 - z3, 3), mp.nstr(c4 - c4claim, 3), mp.nstr(c5 - c5claim, 3)]
