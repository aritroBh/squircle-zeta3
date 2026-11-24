"""Closed forms of the Taylor coefficients a_j(theta) of log R_theta(eps) (see result.md, Theorem 3)."""
import mpmath as mp


def S12(z):
    """Nielsen S_{1,2}(z) = 1/2 int_0^z log^2(1-t) dt/t  (quadrature)."""
    return mp.quad(lambda t: mp.log(1 - t)**2/t, [0, z])/2


def S12_classical(z):
    return mp.zeta(3) - mp.polylog(3, 1 - z) + mp.log(1 - z)*mp.polylog(2, 1 - z) + mp.log(z)*mp.log(1 - z)**2/2


def S22(z):
    """Nielsen S_{2,2}(z) = -1/2 int_0^1 log t log^2(1-zt) dt/t  (quadrature)."""
    return -mp.quad(lambda t: mp.log(t)*mp.log(1 - z*t)**2/t, [0, mp.mpf(1)/2, 1])/2


def closed_a(th):
    z = mp.sin(th/2)**2
    w = (1 - z)/z
    l = mp.log(1 - z)
    u = w*l
    s12 = S12_classical(z)
    a1 = -1 - u
    a2 = mp.zeta(2) + mp.mpf(1)/2 - mp.polylog(2, z) - u**2/2
    a3 = 2*w*s12 - (1 + u**3)/3
    a4 = mp.zeta(4)/2 - 2*S22(z) + 2*w**2*l*s12 + (1 - u**4)/4
    return [a1, a2, a3, a4]
