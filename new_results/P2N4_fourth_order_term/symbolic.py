# P2N4: exact symbolic checks (sympy). Own derivation from the exact formula of thm:P2.
# E_n = 2z q (H-1) + sum_m (-1)^m eps/(m-eps) Lam_m x^m - c2 x^2,  x=1/n, eps=x/(1-x), q=exp(-eps*lam), lam=L-a1,
# H = exp(a2 e^2 + a3 e^3 + a4 e^4 + a5 e^5).
import sympy as sp
x, L, z, a1, a2, a3, a4, a5 = sp.symbols('x L z a1 a2 a3 a4 a5')
c = 1-2*z
P = [sp.Integer(1), c]
for j in range(1, 8): P.append(sp.expand(((2*j+1)*c*P[j]-j*P[j-1])/(j+1)))
Lam = lambda m: sp.factor(sp.expand((P[m-1]-P[m+1])/(2*m+1)))
print("Lam1 =", Lam(1), "; Lam2 =", Lam(2), "; Lam3 =", Lam(3))
eps = x/(1-x)
lam = L - a1
q = sp.exp(-eps*lam)
H = sp.exp(a2*eps**2 + a3*eps**3 + a4*eps**4 + a5*eps**5)
c2 = 2*z*(a2-1+z)
E = 2*z*q*(H-1) + sum((-1)**m*eps/(m-eps)*Lam(m)*x**m for m in range(1, 7)) - c2*x**2
Es = sp.expand(sp.series(E, x, 0, 6).removeO())
co = [sp.expand(Es.coeff(x, k)) for k in range(6)]
print("coeff x^0..x^2:", [sp.simplify(co[k]) for k in range(3)])
d = 2*z*(a2*(2+a1)+a3) - z*(3-z-2*z**2)
print("x^3 - (d - 2 z a2 L):", sp.simplify(co[3] - (d - 2*z*a2*L)))
# Q4 in the compact form of result.md
p4 = -8*z*(1-z) + sp.Rational(3,4)*Lam(2) - Lam(3)/3
print("p4(z) =", sp.expand(p4), " factored:", sp.factor(p4))
Q4 = z*(a2*((lam-3)**2-3) - 2*a3*(lam-3) + a2**2 + 2*a4) + p4
print("x^4 - Q4(compact):", sp.simplify(sp.expand(co[4] - Q4)))
Q4notes = z*(a2*L**2 - (2*a1*a2+6*a2+2*a3)*L + a1**2*a2 + 6*a1*a2 + 2*a1*a3 + a2**2 + 6*a2 + 6*a3 + 2*a4) \
      + sp.Rational(10,3)*z**4 - sp.Rational(11,3)*z**3 + sp.Rational(15,2)*z**2 - sp.Rational(43,6)*z
print("Q4(compact) - Q4(notes/referee):", sp.simplify(sp.expand(Q4 - Q4notes)))
print("x^4 coeff collected in L:", sp.collect(co[4], L))
print("x^5 coeff, L^3 part:", sp.factor(co[5].coeff(L, 3)), "; L^2 part:", sp.factor(co[5].coeff(L, 2)))
# pieces of the proof: formal expansions of P1..P6 separately add to co[3], co[4]
b4 = a4 + a2**2/2
P1 = 2*z*a2*(q*eps**2 - x**2); P2 = 2*z*a3*q*eps**3; P3 = 2*z*b4*q*eps**4
P5 = -4*z*(1-z)*x**3/(1-2*x)
P6 = sum((-1)**m*x/(m-(m+1)*x)*Lam(m)*x**m for m in range(2, 7))
ser = lambda f: sp.expand(sp.series(f, x, 0, 5).removeO())
tot = ser(P1)+ser(P2)+ser(P3)+ser(P5)+ser(P6)
print("pieces x^3 - co3:", sp.simplify(tot.coeff(x,3)-co[3]), "; pieces x^4 - co4:", sp.simplify(tot.coeff(x,4)-co[4]))
print("P1 x^4 coeff - 2 z a2 (3-3lam+lam^2/2):", sp.simplify(ser(P1).coeff(x,4) - 2*z*a2*(3-3*lam+lam**2/2)))
print("P2 x^4 coeff - 2 z a3 (3-lam):", sp.simplify(ser(P2).coeff(x,4) - 2*z*a3*(3-lam)))
print("check eps/(m-eps) = x/(m-(m+1)x):", sp.simplify(eps/(3-eps) - x/(3-4*x)))
print("check T3 = P5:", sp.simplify((-2*z*(1-z)*eps/(1-eps)*x + 2*z*(1-z)*x**2) - P5))
print("check m=2 remainder identity:", sp.simplify(x/(2-3*x) - (x/2 + 3*x**2/4) - (9*x**3/8)/(1-3*x/2)))
print("check m=3 identity:", sp.simplify(1/(3-4*x) - sp.Rational(1,3) - 4*x/(3*(3-4*x))))
print("check (1-2x)^-1-1-2x = 4x^2/(1-2x):", sp.simplify(1/(1-2*x)-1-2*x-4*x**2/(1-2*x)))
