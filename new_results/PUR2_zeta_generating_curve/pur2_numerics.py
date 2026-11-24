"""PUR2 numerics: perimeter by direct u-quadrature of h vs closed form minus tail; convexity data."""
import mpmath as mp
import sys

mp.mp.dps = 60
PI4 = mp.pi / 4


class Curve:
    """phi_c(v) = e^v cosh(c v) (c real) or e^v cos(k v) (c = i k)."""
    def __init__(self, name, c=None, kappa=None):
        self.name = name
        if kappa is None:
            self.c = mp.mpf(c)
            self.a = (1 + self.c) / 2
            self.b = (1 - self.c) / 2
        else:
            self.k = mp.mpf(kappa)
            self.a = mp.mpc(1, self.k) / 2
            self.b = mp.mpc(1, -self.k) / 2
            self.c = None
        self.V = self.Lam_p(mp.mpf(1))

    def phi(self, v):
        if self.c is not None:
            return mp.exp(v) * mp.cosh(self.c * v)
        return mp.exp(v) * mp.cos(self.k * v)

    def dphi(self, v):  # phi'
        if self.c is not None:
            return mp.exp(v) * (mp.cosh(self.c * v) + self.c * mp.sinh(self.c * v))
        return mp.exp(v) * (mp.cos(self.k * v) - self.k * mp.sin(self.k * v))

    def Phi(self, v):
        if self.c is not None:
            if self.c == 1:
                return (mp.exp(2 * v) - 1) / 4 + v / 2
            if self.c == 0:
                return mp.exp(v) - 1
            return ((mp.exp(2 * self.a * v) - 1) / (2 * self.a) + (mp.exp(2 * self.b * v) - 1) / (2 * self.b)) / 2
        k = self.k
        return (mp.exp(v) * (mp.cos(k * v) + k * mp.sin(k * v)) - 1) / (1 + k * k)

    def Lam_p(self, p):
        """inverse of Phi: v with Phi(v) = p."""
        if p == 0:
            return mp.mpf(0)
        if self.c is not None and self.c == 1:     # Lambert W: e^{2v} + 2v = 1 + 4p
            q = 1 + 4 * p
            return (q - mp.lambertw(mp.e**q).real) / 2
        if self.c is not None and self.c == mp.mpf(1) / 2:   # z = e^{v/2}: z^3 + 3z - (4+3p) = 0
            q = 4 + 3 * p
            r = mp.sqrt(q * q / 4 + 1)
            z = mp.cbrt(q / 2 + r) - mp.cbrt(r - q / 2)
            return 2 * mp.log(z)
        v = mp.log(1 + p)
        for _ in range(200):
            dv = (self.Phi(v) - p) / self.phi(v)
            v -= dv
            if abs(dv) < mp.mpf(10)**(-mp.mp.dps + 3):
                break
        return v

    def Lam(self, u):
        return self.Lam_p(1 + mp.sin(u) - mp.cos(u))

    def h(self, u, tau):  # 0 <= u <= pi/4
        return (mp.cos(u) + mp.sin(u)) * (1 + mp.exp(-2 * self.Lam(u) / tau))**tau

    def closed(self, tau):
        a, b = self.a, self.b
        if self.c is not None and self.c == 1:
            main = 8 - 2 * tau * (mp.euler + mp.digamma(1 - tau))
        else:
            R = mp.gamma(1 - a * tau) * mp.gamma(1 - b * tau) * mp.rgamma(1 - tau)
            main = 8 + 2 / (a * b) * (1 - R)
        return main, self.tail(tau)

    def tail(self, tau):  # 8T = 2 tau sum C(tau,n) e^{-2nV/tau} [e^{2aV}/(n-a tau) + e^{2bV}/(n-b tau)]
        V, a, b = self.V, self.a, self.b
        s = 0
        for n in range(1, 2000):
            term = mp.binomial(tau, n) * mp.exp(-2 * n * V / tau) * (mp.exp(2 * a * V) / (n - a * tau) + mp.exp(2 * b * V) / (n - b * tau))
            s += term
            if abs(term) < mp.mpf(10)**(-mp.mp.dps - 5):
                break
        return 2 * tau * s

    def perim_quad(self, tau):
        pts = sorted(set([mp.mpf(0)] + [PI4 * mp.mpf(2)**(-k) for k in range(30, -1, -1)]))
        return 8 * mp.quad(lambda u: self.h(u, tau), pts)

    # ---- convexity ----
    def hpp(self, u, tau):
        """h + h'' on (0, pi/4) from the formula 2 f sigma (M/phi)[-3M' + (M^2/phi)(2K + phi'/phi)]."""
        M, Mp = mp.cos(u) + mp.sin(u), mp.cos(u) - mp.sin(u)
        L = self.Lam(u)
        ph, dph = self.phi(L), self.dphi(L)
        y = 2 * L / tau
        sig = 1 / (1 + mp.exp(y))
        f = (1 + mp.exp(-y))**tau
        K = sig + (1 - sig) / tau
        return 2 * f * sig * (M / ph) * (-3 * Mp + (M * M / ph) * (2 * K + dph / ph))

    def D(self, u):
        """tau-free part: D = -3M' + (M^2/phi)(2 + phi'/phi); returns D*phi/M^2 too."""
        M, Mp = mp.cos(u) + mp.sin(u), mp.cos(u) - mp.sin(u)
        L = self.Lam(u)
        ph, dph = self.phi(L), self.dphi(L)
        Dv = -3 * Mp + (M * M / ph) * (2 + dph / ph)
        return Dv, Dv * ph / (M * M)

    def jump(self, tau):
        V = self.V
        f = (1 + mp.exp(-2 * V / tau))**tau
        sig = 1 / (1 + mp.exp(2 * V / tau))
        return 8 * f * sig / self.phi(V)


def out(*a):
    print(*a)
    sys.stdout.flush()


curves = [Curve('c=1', c=1), Curve('c=1/2', c=mp.mpf(1) / 2), Curve('kappa^2=9/7 (c=i kappa)', kappa=mp.sqrt(mp.mpf(9) / 7)),
          Curve('c=0 (PUR1)', c=0)]
out("V values:")
for C in curves:
    out("  %-26s V = %s   Phi(V)-1 = %s   phi(V) = %s" % (C.name, mp.nstr(C.V, 30), mp.nstr(C.Phi(C.V) - 1, 3), mp.nstr(C.phi(C.V), 15)))
out("  c=1 check V = (5 - W(e^5))/2 :", mp.nstr((5 - mp.lambertw(mp.e**5).real) / 2, 30))

out("\n=== (F) perimeter: direct u-quadrature of h vs closed form minus tail (dps 60) ===")
for C in curves[:3]:
    for tau in ['0.5', '0.2', '0.1', '0.05']:
        tau = mp.mpf(tau)
        I = C.perim_quad(tau)
        main, tl = C.closed(tau)
        diff = main - tl - I
        out("%-26s tau=%-5s int h = %s | closed - tail - int h = %s | tail = %s" % (
            C.name, mp.nstr(tau, 3), mp.nstr(I, 45), mp.nstr(diff, 3), mp.nstr(mp.re(tl), 8)))
# tau = 1
C1 = curves[0]
I = C1.perim_quad(mp.mpf(1))
direct = 8 + 4 * C1.V + 2 * (1 - mp.exp(-2 * C1.V))
out("c=1 tau=1: int h = %s ; 8 + 4V + 2(1-e^{-2V}) - int h = %s" % (mp.nstr(I, 45), mp.nstr(direct - I, 3)))
t1 = 1 - mp.mpf(10)**-25
m_, tl_ = C1.closed(t1)
out("c=1 tau=1-1e-25: closed - tail = %s (limit value above, differs by O(1e-25))" % mp.nstr(m_ - tl_, 30))
for C in curves[1:3]:
    I = C.perim_quad(mp.mpf(1))
    main, tl = C.closed(mp.mpf(1))
    out("%-26s tau=1: int h = %s | closed - tail - int h = %s" % (C.name, mp.nstr(I, 45), mp.nstr(main - tl - I, 3)))

out("\n=== (E) convexity data ===")
for C in curves:
    # formula vs numerical second derivative
    tau = mp.mpf('0.3')
    errs = []
    for u in [mp.mpf('0.05'), mp.mpf('0.4'), mp.mpf('0.75')]:
        num = C.h(u, tau) + mp.diff(lambda s: C.h(s, tau), u, 2)
        errs.append(abs(num - C.hpp(u, tau)) / abs(num))
    hp0 = mp.diff(lambda s: C.h(s, tau), mp.mpf(0), direction=1)
    hpl = mp.diff(lambda s: C.h(s, tau), PI4, direction=-1)
    out("%-26s formula check rel.err %s ; h'(0+) = %s ; jump at pi/4: -2h'(pi/4-) = %s vs 8 f sigma/phi(V) = %s" % (
        C.name, mp.nstr(max(errs), 3), mp.nstr(hp0, 3), mp.nstr(-2 * hpl, 12), mp.nstr(C.jump(tau), 12)))
    Ng = 400
    grid = [PI4 * k / Ng for k in range(1, Ng + 1)]
    Ds = [C.D(u) for u in grid]
    kmin = min(range(Ng), key=lambda k: Ds[k][0])
    out("    min over grid of D = %s at u = %s ; min D phi/M^2 = %s" % (
        mp.nstr(Ds[kmin][0], 8), mp.nstr(grid[kmin], 5), mp.nstr(min(d[1] for d in Ds), 8)))
    for tau in ['1', '0.5', '0.2', '0.1', '0.05']:
        tau = mp.mpf(tau)
        vals = [C.hpp(u, tau) for u in grid]
        k = min(range(Ng), key=lambda i: vals[i])
        out("    tau=%-4s min (h+h'') on grid = %s at u=%s ; jump at pi/4 = %s" % (
            mp.nstr(tau, 3), mp.nstr(vals[k], 8), mp.nstr(grid[k], 4), mp.nstr(C.jump(tau), 8)))
    out("    tau=1.2: (h+h'') at u=1e-6: %s  (expected sign of 1/tau-1 < 0)" % mp.nstr(C.hpp(mp.mpf(10)**-6, mp.mpf('1.2')), 6))
