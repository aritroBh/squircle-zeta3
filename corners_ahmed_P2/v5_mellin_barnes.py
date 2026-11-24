"""V5. Mellin transform of h_c itself, and the Mellin-Barnes form of the exact formula.

  (a) For -1 < Re w < 0:  int_0^oo h_c(s) s^(w-1) ds = -(pi/sin(pi w)) Lambda_w(theta),
      Lambda_nu(theta) := int_{cos th}^1 P_nu(t) dt = 2 sin^2(th/2) 2F1(-nu,1+nu;2;sin^2(th/2)).
  (b) Lambda_nu(theta) = int_{cos th}^1 P_nu(t) dt  (quadrature of the Legendre function) vs the 2F1 form.
  (c) Mellin-Barnes:  L_theta(n) = 2 + (1/2 pi i) int_{(sigma)} (pi/sin pi w) Lambda_w n^-w eps/(eps-w) dw,
      -1 < sigma < 0; checked against the exact formula (V2) at sigma = -1/2.
"""
import mpmath as mp
from common import *

setdps(45)


def Lam(nu, th):
    z = mp.sin(th/2)**2
    return 2*z*mp.hyp2f1(-nu, 1 + nu, 2, z)


angles = [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3), ('3pi/4', 3*mp.pi/4)]

print("V5a: Mellin transform of h_c in the strip -1 < w < 0 (quadrature vs -(pi/sin pi w) Lambda_w)")
for name, th in angles:
    c = mp.cos(th)
    for w in [mp.mpf('-0.25'), mp.mpf('-0.5'), mp.mpf('-0.9')]:
        # int_0^1 h s^(w-1) ds with s = t^(1/(1+w)) ; int_1^oo h s^(w-1) ds with s = t^(-1/(-w))
        a = 1 + w
        b = -w
        f1 = lambda t: (hc(t**(1/a), c)/t**(1/a)) / a if t > 0 else (1 - c)/a
        f2 = lambda t: hc(t**(-1/b), c) / b
        q = mp.quad(f1, [0, mp.mpf(1)/4, 1]) + mp.quad(f2, [0, mp.mpf(1)/4, 1])
        cf = -(mp.pi/mp.sin(mp.pi*w))*Lam(w, th)
        print(f"  theta={name:>6} w={mp.nstr(w,3):>6}  quad={mp.nstr(q,40):>46}  |quad - closed|={mp.nstr(abs(q-cf),3)}")

print("\nV5b: Lambda_nu = int_{cos th}^1 P_nu(t) dt (quadrature of mpmath legendre) vs 2 sin^2(th/2) 2F1(-nu,1+nu;2;sin^2(th/2))")
for name, th in angles:
    for nu in [mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('4')]:
        q = mp.quad(lambda t: mp.legendre(nu, t), [mp.cos(th), 1])
        print(f"  theta={name:>6} nu={mp.nstr(nu,2):>4}  quad={mp.nstr(q,40):>46}  |diff|={mp.nstr(abs(q-Lam(nu, th)),3)}")

print("\nV5c: Mellin-Barnes integral on Re w = -1/2 vs exact formula")
for name, th in angles:
    for n in [mp.mpf(3), mp.mpf(10)]:
        eps = 1/(n - 1)
        sig = mp.mpf(-1)/2

        def f(t):
            w = mp.mpc(sig, t)
            return (mp.pi/mp.sin(mp.pi*w))*Lam(w, th)*mp.power(n, -w)*eps/(eps - w)
        T = 115/(mp.pi - th)
        pts = [0] + [T*k/16 for k in range(1, 17)]
        I = 2*mp.re(mp.quad(f, pts))/(2*mp.pi)   # integrand(-t) = conj(integrand(t))
        mb = 2 + I
        ex = L_exact(n, th)
        print(f"  theta={name:>6} n={mp.nstr(n,3):>4}  MB={mp.nstr(mb,40)}  exact={mp.nstr(ex,40)}  |diff|={mp.nstr(abs(mb-ex),3)}")
