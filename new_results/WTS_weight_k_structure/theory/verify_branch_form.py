"""Numerical (40 digits) check of the algebra used in Sec. 2.4(iii)-(iv):
  w=(1-t)/t, s=l2-l1=log w, sigma=log(i w) (so s = sigma - i pi/2 on the branch with w(i)=-i), D=l1(i)-l2(i)=i pi/2, c=1-eps:
  (a) 1 + w^{-2} e^{2 eps s}  ==  1 - exp(-2 (c*sigma + eps*D))
  (b) t = 1/(1 - i e^sigma),   l1_sigma = d l1/d sigma = i e^sigma/(1 - i e^sigma)   (finite-difference derivative)
at random complex v near i, eps = 0.3 (precision set after imports)."""
import mpmath as mp
mp.mp.dps = 40
print("mp.mp.dps =", mp.mp.dps)
eps = mp.mpf(3)/10; c = 1 - eps
for v in (mp.mpc(0.15, 0.97), mp.mpc(-0.1, 1.04), mp.mpc(0.0, 0.9)):
    t = (v*v + 2*v - 1)/(4*v); w = (1-t)/t
    sigma = mp.log(1j*w)
    s = mp.log(w)                     # principal log; same branch as sigma - i pi/2 near v=i (w near -i)
    lhs = 1 + w**(-2)*mp.e**(2*eps*s)
    D = 1j*mp.pi/2
    rhs = 1 - mp.e**(-2*(c*sigma + eps*D))
    ta = 1/(1 - 1j*mp.e**sigma)
    h = mp.mpf(10)**(-15)
    l1 = lambda sg: -mp.log(1 - 1j*mp.e**sg)
    dl1 = (l1(sigma+h) - l1(sigma-h))/(2*h)
    print("v=", v, " |(a) diff| =", mp.nstr(abs(lhs-rhs), 3), " |s-(sigma-i pi/2)| =", mp.nstr(abs(s-(sigma-1j*mp.pi/2)), 3),
          " |(b) t diff| =", mp.nstr(abs(t-ta), 3), " |(b) l1_sigma diff| =", mp.nstr(abs(dl1 - 1j*mp.e**sigma/(1-1j*mp.e**sigma)), 3))
