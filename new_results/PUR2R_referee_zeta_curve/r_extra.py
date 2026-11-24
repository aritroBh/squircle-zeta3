"""Referee PUR2R: Remark A' (continuation past the poles) for real c = 2, where Theorem A's
direct range is tau < 2/3 and Gamma(1 - a tau) has a pole at tau = 2/3; and kappa^2 = 9/7 at tau = 1.5."""
import mpmath as mp
from r_curve import Body
mp.mp.dps = 50
c = mp.mpf(2); a = (1 + c) / 2; b = (1 - c) / 2
phi = lambda v: mp.exp(v) * mp.cosh(c * v)
Phi = lambda v: ((mp.exp(2 * a * v) - 1) / (2 * a) + (mp.exp(2 * b * v) - 1) / (2 * b)) / 2
def inv(p):
    return mp.findroot(lambda v: Phi(v) - p, (mp.mpf(0), mp.mpf(2)), solver='anderson') if p > 0 else mp.mpf(0)
V = inv(mp.mpf(1))
def perim(tau):
    return 8 * mp.quad(lambda u: (mp.cos(u) + mp.sin(u)) * (1 + mp.exp(-2 * inv(1 + mp.sin(u) - mp.cos(u)) / tau)) ** tau, [0, 0.1, 0.4, mp.pi / 4])
def closed(tau):
    R = mp.gamma(1 - a * tau) * mp.gamma(1 - b * tau) * mp.rgamma(1 - tau)
    main = 8 + (2 / (a * b)) * (1 - R)
    E = 0
    for n in range(1, 400):
        E += 2 * tau * mp.binomial(tau, n) * mp.exp(-2 * n * V / tau) * (mp.exp(2 * a * V) / (n - a * tau) + mp.exp(2 * b * V) / (n - b * tau))
    return main - E
mp.mp.dps = 30
for tau in [mp.mpf('0.5'), mp.mpf('0.8'), mp.mpf(2) / 3 + mp.mpf('1e-12'), mp.mpf('1')]:
    P = perim(tau); C = closed(tau)
    print("c=2 tau=%s: quad = %s ; closed - E = %s ; diff = %s" % (mp.nstr(tau, 14), mp.nstr(P, 25), mp.nstr(C, 25), mp.nstr(P - C, 3)))
B = Body('imag', mp.mpf(9) / 7)
tau = mp.mpf('1.5')
P = 8 * mp.quad(lambda u: B.h_plus(u, tau), [0, 0.1, 0.4, mp.pi / 4])
print("kappa^2=9/7 tau=1.5 (signed length): quad - (closed - E) =", mp.nstr(P - mp.re(B.closed_minus_tail(tau)[0]), 3))
