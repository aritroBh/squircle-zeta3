"""RULQR script 5: independent evaluation of the truncation bound of RULQ Theorem B2,
|R_K(p)| <= R0 + R1 + R2, in two ways on the same x- and y-grids:
  'mid'  : midpoint Riemann sums (ordinary numerics, to compare with RULQ's ru4 numbers),
  'cert' : rigorous upper Riemann sums in Arb ball arithmetic (every factor bounded above on each cell:
           Gamma(N+1, z) and (x^2+r^2)^{-N/2} are decreasing, so they are evaluated at the cell end that maximises them;
           m1 is evaluated on the whole cell as a ball), plus closed-form tail bounds (REFEREE.md Sec. 5).
Usage: rr5_bound.py p K rho delta rho2 delta2 [grid factor c]"""
import sys
from flint import arb, acb, ctx
ctx.dps = 30
p, K = arb(sys.argv[1]), int(sys.argv[2])
rho, delta, rho2, delta2 = (arb(s) for s in sys.argv[3:7])
cfac = float(sys.argv[7]) if len(sys.argv) > 7 else 0.02
h = 1 / p; pi = arb.pi(); hf = float(h.mid())
assert (h < rho) and ((1 + rho) * h < 1) and (h < rho2) and (rho2 < 1) and (delta <= pi / 2) and (delta2 <= pi), 'parameters'
fact = [arb(1)]
for i in range(1, K + 3): fact.append(fact[-1] * i)

def gam_list(z, Nmax):
    """[h^{N+1} Gamma(N+1, z) for N = 0..Nmax], z >= 0 an arb (finite-sum formula, exact for integer N)."""
    ez = (-z).exp(); out = []; s = arb(0); term = arb(1)
    for N in range(Nmax + 1):
        if N > 0: term = term * z / N
        s = s + term
        out.append(fact[N] * h ** (N + 1) * ez * s)
    return out

def make_grid(lo, hi, feats, small, c):
    """points from lo to hi; cell width c*max(small, distance to nearest feature), capped by c*hf beyond the first feature
    that is followed by Gamma decay (handled by caller through 'small')."""
    pts = [lo]; x = lo
    while x < hi:
        d = min(abs(x - f) for f in feats)
        w = c * max(small(x), d)
        x = min(x + w, hi); pts.append(x)
    return pts

def xpart(r, M, a, Ns, ws, mfun, cert):
    """sum_j ws[j] * (1/(pi r)) int_R m(x) (x^2+r^2)^{-N_j/2} h^{N_j+1} Gamma(N_j+1,(x-r)^+/h) dx, plus tails."""
    rf = float(r.mid()); Nmax = max(Ns)
    XL = 50 * rf; XR = rf + hf * (3 * K + 80)
    small = lambda x: float(delta.mid()) if x < rf else hf
    pts = make_grid(-XL, XR, [0.0, rf], small, cfac)
    tot = arb(0)
    for xa, xb in zip(pts[:-1], pts[1:]):
        A, B = arb(xa), arb(xb); w = B - A
        if cert:
            xball = arb((xa + xb) / 2, (xb - xa) / 2)
            mval = arb(mfun(xball).upper())
            xmin = arb(0) if xa <= 0 <= xb else (B if xb < 0 else A)
            zlo = A - r; zlo = arb(0) if zlo.upper() <= 0 else arb(zlo.lower()) / h
        else:
            xm = (A + B) / 2; mval = mfun(xm); xmin = xm
            zlo = xm - r; zlo = arb(0) if zlo.mid() <= 0 else zlo / h
        g = gam_list(zlo, Nmax); base = xmin * xmin + r * r
        s = sum((ws[i] * base ** (arb(-Ns[i]) / 2) * g[Ns[i]] for i in range(len(Ns))), arb(0))
        val = mval * s * w
        tot += arb(val.upper()) if cert else val
    # tails: x < -XL (m <= M, Gamma <= N!, (x^2+r^2)^{-N/2} <= |x|^{-N}); x > XR (m <= M e^{a x}, theta = 9/10)
    th = arb(9) / 10; b = th / h - a
    tl = sum((ws[i] * M * fact[N] * h ** (N + 1) * arb(XL) ** (1 - N) / (N - 1) for i, N in enumerate(Ns)), arb(0))
    tr = sum((ws[i] * M * fact[N] * h ** (N + 1) * (1 - th) ** (-N - 1) * arb(XR) ** (-N) for i, N in enumerate(Ns)), arb(0)) \
         * (th * r / h - b * XR).exp() / b
    return (tot + tl + tr) / (pi * r)

def ypart(r, M, a, Ns, ws, cert):
    """sum_j ws[j] * (M/pi) e^{a r} int_0^oo e^{-y(1/h-a)} y^N (y+r)^{-N} dy, upper sums + tail."""
    b = 1 / h - a; Y = float(((arb(60) + 2 * K) / b).mid()) + 1
    pts = make_grid(0.0, Y, [0.0], lambda y: hf, cfac)
    tot = arb(0)
    for ya, yb in zip(pts[:-1], pts[1:]):
        A, B = arb(ya), arb(yb)
        if cert: e = (-b * A).exp(); q = B / (B + r)
        else: m = (A + B) / 2; e = (-b * m).exp(); q = m / (m + r)
        s = sum((ws[i] * q ** Ns[i] for i in range(len(Ns))), arb(0))
        val = e * s * (B - A)
        tot += arb(val.upper()) if cert else val
    tot += sum(ws, arb(0)) * (-b * Y).exp() / b
    return M / pi * (a * r).exp() * tot

def bound(cert):
    R0 = 2 * (h * arb(2).log()) ** (K + 1) / fact[K + 1]
    # n = 1 part
    r1 = pi - delta; s1 = 2 ** arb(0.5) * (delta / 2).sin(); Lam = pi + (2 / s1).log()
    M1 = (rho * Lam).exp() / s1; a1 = 1 + rho
    def m1(x):
        z = 1 + acb(x, r1).exp()
        return (2 * x).exp() * (rho * abs(z.log())).exp() / abs(z)
    Ns1 = [K - 2 - j for j in range(0, K - 3)]; ws1 = [(h / rho) ** j for j in range(0, K - 3)]
    R1 = xpart(r1, M1, a1, Ns1, ws1, m1, cert) + ypart(r1, M1, a1, Ns1, ws1, cert)
    R1 = h ** 2 * (R1 + 2 * M1 * h * (h / rho) ** (K - 3) / ((1 - h / rho) * (1 - a1 * h)))
    # n >= 2 part
    r2 = 2 * pi - delta2; s2 = 2 ** arb(0.5) * (delta2 / 4).sin(); a2 = 1 + rho2 / 2
    M2 = rho2 ** 2 * (rho2 * (1 + pi)).exp() * s2 ** (-1 - rho2) * ((2 - rho2).zeta() - 1)
    m2 = lambda x: M2 * (a2 * (x if x.lower() > 0 else (arb(0) if x.upper() < 0 else arb(0).union(x)))).exp() if cert else \
         M2 * (a2 * (x if x.mid() > 0 else arb(0))).exp()
    Ns2 = [K - j for j in range(2, K - 1)]; ws2 = [(h / rho2) ** j for j in range(2, K - 1)]
    R2 = xpart(r2, M2, a2, Ns2, ws2, m2, cert) + ypart(r2, M2, a2, Ns2, ws2, cert)
    R2 = R2 + 2 * M2 * h * (h / rho2) ** (K - 1) / ((1 - h / rho2) * (1 - a2 * h))
    return R0 + R1 + R2, R0, R1, R2

for mode in ('mid', 'cert'):
    tot, R0, R1, R2 = bound(mode == 'cert')
    f = lambda v: arb(v.upper()).str(5, radius=False) if mode == 'cert' else v.str(5, radius=False)
    print('p=%s K=%d rho=%s delta=%s rho2=%s delta2=%s c=%g  %-4s  bound %s   (R0 %s, R1 %s, R2 %s)' % (
        sys.argv[1], K, *sys.argv[3:7], cfac, mode, f(tot), f(R0), f(R1), f(R2)), flush=True)
