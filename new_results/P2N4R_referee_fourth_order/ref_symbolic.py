# P2N4R referee: own sympy derivation of the n^-3, n^-4, n^-5 coefficients of E_n from the exact formula of thm:P2.
# Built independently of P2N4/symbolic.py: Legendre polynomials from sympy.legendre, x=1/n, log n = L symbol,
# g(eps) = exp(sum_j a_j eps^j) is the DEFINITION of the a_j (log of R_theta(eps)/(2z)), n^-eps = exp(-eps L).
import sympy as sp
x, L, z = sp.symbols('x L z', positive=True)
a = sp.symbols('a1:7')
a1, a2, a3, a4, a5, a6 = a
c = 1 - 2*z
Lam = lambda m: sp.expand((sp.legendre(m-1, c) - sp.legendre(m+1, c))/(2*m+1))
eps = x/(1-x)
N = 7  # expand to x^6 inclusive
logg = sum(a[j-1]*eps**j for j in range(1, 7))
# full exact formula, then subtract the displayed part of thm:P2asym
Lfull = 1 + c + 2*z*sp.exp(logg - eps*L) + sum((-1)**m*eps/(m-eps)*Lam(m)*x**m for m in range(1, N+1))
c2 = 2*z*(a2 - 1 + z)
shown = 1 + c + 2*z*sp.exp(a1*eps - eps*L) + c2*x**2
E = Lfull - shown
Es = sp.expand(sp.series(E, x, 0, 6).removeO())
co = {k: sp.expand(Es.coeff(x, k)) for k in range(6)}
out = []
def p(*s):
    print(*s); out.append(' '.join(map(str, s)))
p("coeff x^0, x^1, x^2:", co[0], co[1], co[2])
d = 2*z*(a2*(2+a1)+a3) - z*(3-z-2*z**2)
p("x^3 coeff - (d - 2 z a2 L) =", sp.simplify(co[3] - (d - 2*z*a2*L)))
lam = L - a1
p4 = -sp.Rational(1, 6)*z*(1-z)*(20*z**2-2*z+43)
Q4 = z*(a2*((lam-3)**2-3) - 2*a3*(lam-3) + a2**2 + 2*a4) + p4
p("x^4 coeff - Q4(claimed compact form) =", sp.simplify(co[4] - Q4))
Q4exp = z*a2*L**2 - 2*z*(a1*a2+3*a2+a3)*L + z*(a1**2*a2+6*a1*a2+6*a2+2*a1*a3+6*a3+a2**2+2*a4) \
        + sp.Rational(10,3)*z**4 - sp.Rational(11,3)*z**3 + sp.Rational(15,2)*z**2 - sp.Rational(43,6)*z
p("x^4 coeff - Q4(claimed expanded form) =", sp.simplify(co[4] - Q4exp))
p("x^4 coeff does not involve a5, a6:", all(sp.diff(co[4], s) == 0 for s in (a5, a6)))
Q5 = sp.collect(co[5], L)
p("x^5 coeff, L^3 part:", sp.factor(co[5].coeff(L, 3)))
p("x^5 coeff, L^2 part:", sp.factor(co[5].coeff(L, 2)))
p("x^5 coeff, L^1 part:", sp.factor(co[5].coeff(L, 1)))
p("x^5 coeff, L^0 part:", sp.expand(co[5].coeff(L, 0)))
p("Lambda_2 =", sp.factor(Lam(2)), "; Lambda_3 =", sp.factor(Lam(3)))
p("p4 check (-8z(1-z)+3/4 Lam2-1/3 Lam3) - p4 =", sp.simplify(-8*z*(1-z) + sp.Rational(3,4)*Lam(2) - Lam(3)/3 - p4))
# Step-6 identity of the proof
p("2a2(3-3lam+lam^2/2) - a2((lam-3)^2-3) =", sp.simplify(2*a2*(3-3*lam+lam**2/2) - a2*((lam-3)**2-3)))
# b4 = coefficient of eps^4 of H = exp(sum_{j>=2} a_j eps^j)
e = sp.symbols('e')
Hs = sp.expand(sp.series(sp.exp(a2*e**2+a3*e**3+a4*e**4+a5*e**5), e, 0, 6).removeO())
p("H coefficients e^0..e^5:", [Hs.coeff(e, k) for k in range(6)])
# Step 0(i): r2,r3,r4 Taylor coefficients and Lagrange constants
for k, J in ((2, 2), (3, 1), (4, 0), (5, 0)):
    p("(1-x)^-%d: binom(k+J,J+1) =" % k, sp.binomial(k+J, J+1), "; series:", sp.series((1-x)**(-k), x, 0, J+2))
# The pieces P1..P6 of the proof, formal expansion to x^5, compared with the full coefficients
q = sp.exp(-eps*lam)
b4 = a4 + a2**2/2
P1 = 2*z*a2*(q*eps**2 - x**2); P2 = 2*z*a3*q*eps**3; P3 = 2*z*b4*q*eps**4
P4 = 2*z*q*(a5 + a2*a3)*eps**5   # formal leading part of rho_5
P5 = -4*z*(1-z)*x**3/(1-2*x)
P6 = sum((-1)**m*x/(m-(m+1)*x)*Lam(m)*x**m for m in range(2, N+1))
tot = sp.expand(sum(sp.series(P, x, 0, 6).removeO() for P in (P1, P2, P3, P4, P5, P6)))
for k in (3, 4, 5):
    p("pieces x^%d - full x^%d =" % (k, k), sp.simplify(tot.coeff(x, k) - co[k]))
# P1's x^5 coefficient: its L^3 part is the only L^3 in co[5]
P1s = sp.expand(sp.series(P1, x, 0, 6).removeO())
p("P1 x^5 coefficient in L (lam = L - a1):", sp.collect(sp.expand(P1s.coeff(x, 5)), L))
p("P5 identity T3:", sp.simplify((-2*z*(1-z)*eps/(1-eps)*x + 2*z*(1-z)*x**2) - P5))
open('ref_symbolic.log', 'w').write('\n'.join(out) + '\n')
