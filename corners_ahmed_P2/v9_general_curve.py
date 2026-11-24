"""V9. General monotone corner curve x u + f(x) v (f increasing, f(0)=0, f(1)=1):
   2 - L = int_0^1 h_c(f'(x)) dx = (1/2 pi i) int_{(sigma)} H_theta(w) Mf(w) dw,   -1 < sigma < 0,
   H_theta(w) = -(pi/sin pi w) Lambda_w(theta),   Mf(w) = int_0^1 f'(x)^(-w) dx.
Test family: f(x) = (e^{a x} - 1)/(e^a - 1), Mf(w) = (a/(e^a-1))^(-w) (1 - e^{-a w})/(a w)."""
import mpmath as mp
from common import *
setdps(40)
for name, th in [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3), ('3pi/4', 3*mp.pi/4)]:
    c = mp.cos(th); z = mp.sin(th/2)**2
    for a in [mp.mpf(3), mp.mpf(12)]:
        A = a/(mp.e**a - 1)
        direct = mp.quad(lambda x: hc(A*mp.e**(a*x), c), mp.linspace(0, 1, 9))
        def f(t):
            w = mp.mpc(-0.5, t)
            H = -(mp.pi/mp.sin(mp.pi*w))*2*z*mp.hyp2f1(-w, 1 + w, 2, z)
            M = A**(-w)*(1 - mp.e**(-a*w))/(a*w)
            return H*M
        T = 110/(mp.pi - th)
        mb = 2*mp.re(mp.quad(f, [0] + [T*k/16 for k in range(1, 17)]))/(2*mp.pi)
        print(f"  theta={name:>6} a={mp.nstr(a,3):>4}: 2-L direct={mp.nstr(direct,35)}  MB={mp.nstr(mb,35)}  |diff|={mp.nstr(abs(direct-mb),3)}")
