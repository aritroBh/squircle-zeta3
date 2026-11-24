# residue_k2.py -- Laurent residues at sigma=0 of K_n(sigma) b^(-sigma),
#   K_n(sigma) = -(1/2) Gamma(-sigma/2) Gamma^{(n)}(sigma) (1-sigma)^n,  n = 0,1,2,
# using Gamma(x) = Gamma(1+x)/x and log Gamma(1+x) = -gamma x + sum_k (-1)^k zeta(k) x^k / k.
import sympy as sp
s, l = sp.symbols('sigma ell')
g = sp.EulerGamma
N = 7
lg = -g*s + sum((-1)**k*sp.zeta(k)*s**k/k for k in range(2, N))
G1 = sp.series(sp.exp(lg), s, 0, N).removeO()          # Gamma(1+s)
Gs = sp.expand(G1/s)                                     # Gamma(s) Laurent (truncated)
Ghalf = sp.expand(G1.subs(s, -s/2)/(-s/2))               # Gamma(-s/2)
for n in range(3):
    Gn = sp.diff(Gs, s, n)
    expr = sp.expand(-sp.Rational(1,2)*Ghalf*Gn*(1-s)**n*sp.series(sp.exp(-s*l), s, 0, N).removeO())
    res = sp.expand(expr.coeff(s, -1))
    print("n=%d: Res =" % n, sp.simplify(res))
