"""Referee PUR2R: own implementation of the body K_{tau,c} (c = 1 and c = i*kappa)."""
import mpmath as mp

class Body:
    def __init__(self, kind, kappa2=None):
        self.kind = kind                       # 'c1' or 'imag'
        if kind == 'imag':
            self.k = mp.sqrt(mp.mpf(kappa2))
        self.V = self.Lam_of_p(mp.mpf(1))
        self.cache = {}

    def phi(self, v):
        return mp.exp(v) * (mp.cosh(v) if self.kind == 'c1' else mp.cos(self.k * v))

    def dphi(self, v):
        if self.kind == 'c1':
            return mp.exp(v) * (mp.cosh(v) + mp.sinh(v))
        return mp.exp(v) * (mp.cos(self.k * v) - self.k * mp.sin(self.k * v))

    def Phi(self, v):   # int_0^v phi, by antiderivative
        if self.kind == 'c1':
            return (mp.expm1(2 * v)) / 4 + v / 2
        k = self.k
        return (mp.exp(v) * (mp.cos(k * v) + k * mp.sin(k * v)) - 1) / (1 + k * k)

    def Lam_of_p(self, p):
        """Newton with bisection safeguard on [0, vmax], Phi increasing there."""
        if p == 0:
            return mp.mpf(0)
        lo = mp.mpf(0)
        hi = mp.mpf(2) if self.kind == 'c1' else mp.pi / (2 * self.k)
        v = mp.log(1 + p) / 1.3
        tol = mp.mpf(10) ** (-mp.mp.dps + 4)
        for _ in range(500):
            F = self.Phi(v) - p
            if F > 0: hi = v
            else: lo = v
            vn = v - F / self.phi(v)
            if not (lo < vn < hi):
                vn = (lo + hi) / 2
            if abs(vn - v) < tol:
                v = vn
                break
            v = vn
        return v

    def Lam(self, u):     # u in [0, pi/4]
        key = (u, mp.mp.dps)
        if key not in self.cache:
            self.cache[key] = self.Lam_of_p(1 + mp.sin(u) - mp.cos(u))
        return self.cache[key]

    # one window, 0 <= u <= pi/4
    def h_plus(self, u, tau):
        L = self.Lam(u)
        return (mp.cos(u) + mp.sin(u)) * (1 + mp.exp(-2 * L / tau)) ** tau

    def hp_plus(self, u, tau):   # h' from the right-window formula, via chain rule
        L = self.Lam(u); M = mp.cos(u) + mp.sin(u); Mp = mp.cos(u) - mp.sin(u)
        f = (1 + mp.exp(-2 * L / tau)) ** tau
        sig = 1 / (1 + mp.exp(2 * L / tau))
        Lp = M / self.phi(L)
        return Mp * f - M * 2 * f * sig * Lp

    def hhpp_plus(self, u, tau):  # h + h'' by my own derivation (see REFEREE.md)
        L = self.Lam(u); M = mp.cos(u) + mp.sin(u); Mp = mp.cos(u) - mp.sin(u)
        ph = self.phi(L); dph = self.dphi(L)
        f = (1 + mp.exp(-2 * L / tau)) ** tau
        sig = 1 / (1 + mp.exp(2 * L / tau))
        K = sig + (1 - sig) / tau
        return 2 * f * sig * (M / ph) * (-3 * Mp + (M * M / ph) * (2 * K + dph / ph))

    # whole circle
    def _red(self, u):
        j = mp.nint(u / (mp.pi / 2))
        return u - j * mp.pi / 2

    def h(self, u, tau):
        t = self._red(u)
        return self.h_plus(abs(t), tau)

    def hp(self, u, tau, side=0):
        """h'(u); at a window boundary choose side=-1 (left limit) or +1 (right limit)."""
        t = self._red(u)
        if abs(abs(t) - mp.pi / 4) < mp.mpf(10) ** (-mp.mp.dps + 5) and side != 0:
            # left limit at pi/4 + j pi/2 belongs to the window with t -> +pi/4
            return self.hp_plus(mp.pi / 4, tau) * (1 if side < 0 else -1)
        return self.hp_plus(abs(t), tau) * (1 if t >= 0 else -1)

    # closed forms
    def closed_minus_tail(self, tau):
        tau = mp.mpf(tau); V = self.V
        if self.kind == 'c1':
            main = 8 - 2 * tau * (mp.euler + mp.digamma(1 - tau))
            a, b = mp.mpf(1), mp.mpf(0)
        else:
            a = mp.mpc(1, self.k) / 2; b = mp.mpc(1, -self.k) / 2
            R = mp.gamma(1 - a * tau) * mp.gamma(1 - b * tau) / mp.gamma(1 - tau)
            main = 8 + (2 / (a * b)) * (1 - R)
        E = 0; n = 1
        while True:
            term = 2 * tau * mp.binomial(tau, n) * mp.exp(-2 * n * V / tau) * (
                mp.exp(2 * a * V) / (n - a * tau) + mp.exp(2 * b * V) / (n - b * tau))
            E += term
            if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 5) and n > 5:
                break
            n += 1
        return main - E, E, main
