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


