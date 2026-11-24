# R1: independent symbolic checks of (a) q(p), q^{-p} expansion, (b) order-1/2 forcing of the ODE,
# (c) N(p) Gamma formula derivatives N'(2), N''(2), N'''(2), (d) N(p) vs direct quadrature at p=3.
import sympy as sp
p, e, th = sp.symbols('p e theta', positive=True)
# (a) Phi_p = (pi/2) F_p/F_p(1), F_p(1) = B(1/p,1-1/p)/p = pi/(p sin(pi/p)); S = Phi^{-1}; with sin_p = F_p^{-1}:
# S(theta) = sin_p(2 F_p(1) theta/pi) so  sin_p'(x) = S'(theta) * pi/(2F_p(1)) = q S',  q = p sin(pi/p)/2.
q = p*sp.sin(sp.pi/p)/2
Q = sp.series((q**(-p)).subs(p, 2+e), e, 0, 4).removeO()
print('q(2) =', q.subs(p, 2), '  q^{-p} =', sp.simplify(sp.expand(Q)))
# (b) S'' = -q^{-p} (S/S')^{p-2} S ; expand S = sin + e s1 + e^2/2 s2
s1, s2 = sp.Function('s1')(th), sp.Function('s2')(th)
S = sp.sin(th) + e*s1 + e**2/2*s2
rhs = -Q*sp.exp(e*sp.log(S/sp.diff(S, th)))*S
ser = sp.series(rhs, e, 0, 3).removeO()
c1 = sp.simplify(ser.coeff(e, 1)); c2 = sp.simplify(2*ser.coeff(e, 2))
L0 = sp.log(sp.tan(th))
F1 = (1 - L0)*sp.sin(th)
F2 = -2*L0*s1 + 2*sp.tan(th)*sp.diff(s1, th) - L0**2*sp.sin(th) + 2*L0*sp.sin(th) - (sp.Rational(1, 2) + sp.pi**2/8)*sp.sin(th)
d1 = sp.simplify(sp.expand_log((c1 + s1) - F1, force=True))
d2 = sp.simplify(sp.expand_log(sp.expand((c2 + s2) - F2), force=True))
print('order1: s1''+s1 - (1-L0)sin =', d1)
print('order2: s2''+s2 - F2 =', sp.simplify(sp.expand(sp.expand_log(sp.expand_trig(d2), force=True))))
# (c) N(p)
N = 2*sp.sin(sp.pi/p)/sp.pi*sp.gamma(3/p)*sp.gamma(1-1/p)/sp.gamma(1+2/p)
Ns = sp.series(N.subs(p, 2+e), e, 0, 4).removeO()
for k in range(4):
    v = sp.nsimplify(sp.simplify(sp.factorial(k)*Ns.coeff(e, k)).rewrite(sp.zeta))
    print('N^(%d)(2) =' % k, sp.simplify(sp.expand(sp.factorial(k)*Ns.coeff(e, k))))
N3 = (87*sp.zeta(3) + 16*sp.log(2)**3 - 96*sp.log(2)**2 + 72*sp.log(2) + 10*sp.pi**2*sp.log(2) - 20*sp.pi**2)/16
print("N'''(2) - claimed =", sp.N(6*Ns.coeff(e, 3) - N3, 50))
N2 = sp.log(2)**2 - 3*sp.log(2) + 5*sp.pi**2/24
print("N''(2) - claimed =", sp.N(2*Ns.coeff(e, 2) - N2, 50))
# (d) N(3) by quadrature of (4/pi) int_0^{pi/2} S^2 dtheta = (4/pi) int_0^1 y^2 Phi'(y) dy
from mpmath import mp, quad, hyp2f1, pi, sin, mpf, gamma, nstr
mp.dps = 30
pp = mpf(3); F1v = pi/(pp*sin(pi/pp))
val = 4/pi*quad(lambda y: y**2*(pi/2)/F1v*(1-y**pp)**(-1/pp), [0, 1])
print('N(3) quad - Gamma formula =', nstr(val - 2*sin(pi/pp)/pi*gamma(3/pp)*gamma(1-1/pp)/gamma(1+2/pp), 3))
