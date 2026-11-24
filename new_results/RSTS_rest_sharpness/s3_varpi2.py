"""Symbolic first correction of the q=2 saddle:  Z_k = Zpred*(1 + varpi2/sqrt(N) + O(1/N)).
Z_k = int H E A2_k ds; s = 2 lam w, lam = sqrt(N/pi), Lam = sqrt(pi N) = pi*lam.
Integrand = C * w^{-5/2} e^{pi/(2w^2)} e^{Lam phi0(w)} [1 + (1/lam) c1(w) + O(lam^-2)],
phi0 = i w - 1/w, c1(w) = 1/(8w) - pi/(3 w^3)   (from H, (2s+2)/s^2 and (1+x)^-N to O(x^3)).
Laplace at w0 = e^{i pi/4}: second-order term via Gaussian moments."""
import sympy as sp
t, L = sp.symbols('t Lambda', positive=True)
w0 = sp.exp(sp.I*sp.pi/4)
w = w0 + t
phi = sp.I*w - 1/w
g = w**sp.Rational(-5, 2)*sp.exp(sp.pi/(2*w**2))
c1 = 1/(8*w) - sp.pi/(3*w**3)
ser = lambda e, n: sp.series(e, t, 0, n).removeO()
phis = sp.expand(ser(phi, 5))
a = [sp.simplify(phis.coeff(t, j)) for j in range(5)]
assert sp.simplify(a[1]) == 0
gs = sp.expand(ser(g, 3)); gc = [sp.simplify(gs.coeff(t, j)) for j in range(3)]
g0 = gc[0]
a2, a3, a4 = a[2], a[3], a[4]
# moments: int t^{2m} e^{L a2 t^2} dt = Gamma(m+1/2) (-L a2)^{-m-1/2}; ratio to m=0: mom(m) = Gamma(m+1/2)/Gamma(1/2) (-L a2)^{-m}
mom = lambda m: sp.gamma(m+sp.Rational(1, 2))/sp.sqrt(sp.pi)*(-L*a2)**(-m)
# relative correction from Laplace (order 1/L):
rel_L = (gc[2]*mom(1) + gc[1]*L*a3*mom(2) + g0*(L*a4*mom(2) + L**2*a3**2/2*mom(3)))/g0
rel_L = sp.simplify(sp.expand(rel_L))
# relative correction from c1 at order 1/lam = pi/L:
rel_c = sp.pi/L*c1.subs(t, 0)
tot = sp.simplify(sp.expand((rel_L + rel_c)*L))  # coefficient of 1/Lambda
# 1/Lambda = 1/(sqrt(pi) sqrt(N))  => varpi2 = tot/sqrt(pi)
varpi2 = sp.nsimplify(sp.simplify(tot/sp.sqrt(sp.pi)))
re_, im_ = [sp.simplify(sp.expand_complex(x)) for x in (sp.re(sp.expand_complex(tot/sp.sqrt(sp.pi))), sp.im(sp.expand_complex(tot/sp.sqrt(sp.pi))))]
print("varpi2 = ", sp.simplify(re_), " + i*(", sp.simplify(im_), ")")
print("numeric:", sp.N(re_, 40), sp.N(im_, 40))
