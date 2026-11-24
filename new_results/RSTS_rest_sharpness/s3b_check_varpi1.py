"""Control: the same symbolic machinery applied to the q=1 main integral must reproduce the paper's varpi_1."""
import sympy as sp
t, L = sp.symbols('t Lambda', positive=True)
w0 = sp.exp(sp.I*sp.pi/4); w = w0 + t
phi = sp.I*w - 1/w
g = w**sp.Rational(-5, 2)*sp.exp(sp.pi/(4*w**2))
c1 = -sp.Rational(3, 4)/w - sp.pi/(6*w**3)        # in units of 1/mu, mu = sqrt(2N/pi), 1/mu = pi/(2 Lambda1)
ser = lambda e, n: sp.series(e, t, 0, n).removeO()
phis = sp.expand(ser(phi, 5)); a = [phis.coeff(t, j) for j in range(5)]
gs = sp.expand(ser(g, 3)); gc = [gs.coeff(t, j) for j in range(3)]
a2, a3, a4 = a[2], a[3], a[4]
mom = lambda m: sp.gamma(m+sp.Rational(1, 2))/sp.sqrt(sp.pi)*(-L*a2)**(-m)
rel_L = (gc[2]*mom(1) + gc[1]*L*a3*mom(2) + gc[0]*(L*a4*mom(2) + L**2*a3**2/2*mom(3)))/gc[0]
tot = sp.expand((rel_L + sp.pi/(2*L)*c1.subs(t, 0))*L)   # coefficient of 1/Lambda1, Lambda1 = sqrt(pi N/2)
v = sp.expand_complex(tot*sp.sqrt(2)/sp.sqrt(sp.pi))
varpi1 = ((24+12*sp.pi+sp.pi**2) + sp.I*(24-12*sp.pi+sp.pi**2))/(48*sp.sqrt(sp.pi))
print("derived:", sp.N(v, 30)); print("paper  :", sp.N(varpi1, 30)); print("diff:", sp.N(v-varpi1, 30))
