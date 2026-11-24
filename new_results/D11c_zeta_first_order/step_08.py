import sympy as sp
from sympy import symbols, Rational, factorial, binomial, sqrt, exp, series, simplify

NMAX = 16
t = symbols('t')
Rt = series(sqrt((1+exp(4*t))/2), t, 0, NMAX+1).removeO()
V = {n: sp.simplify(factorial(n)*Rt.coeff(t, n)) for n in range(NMAX+1)}

kappa_exact = {3: Rational(1), 5: Rational(-3,4), 7: Rational(129,16), 9: Rational(-22805,192),
               11: Rational(805173,256), 13: Rational(-131893755,1024), 15: Rational(154940467501,20480)}

# c_k exact
c = {0: Rational(1), 1: Rational(1,4), 2: Rational(-1,96), 3: Rational(19,5760), 4: Rational(-559,645120),
     5: Rational(29161,116121600), 6: Rational(-2368081,30656102400), 7: Rational(276580459,11158821273600)}

print(f"{'n':>4} {'kappa exact':>22} {'V-form':>22} {'c_k-sum form':>22}  match")
for n in range(3, 16, 2):
    m = (n-1)//2
    kV = sp.Rational(4*(-1)**n, n) * (1 - sp.Rational(2)**(1-n) * V[n])
    kc = sp.Rational(4*(-1)**n, n) * (1 - 2*sum(c[k]*factorial(n)/factorial(n-2*k)*sp.Rational(2)**(2*k-n) for k in range(0, m+1)))
    ok = (sp.simplify(kV - kappa_exact[n])==0) and (sp.simplify(kc - kappa_exact[n])==0)
    print(f"{n:>4} {str(kappa_exact[n]):>22} {str(sp.simplify(kV)):>22} {str(sp.simplify(kc)):>22}  {ok}")

print("\nV_n odd values:", {n: V[n] for n in range(1, 16, 2)})
# factor denominators/numerators curiosity
from sympy import factorint
for n in range(3, 16, 2):
    print(n, "kappa =", kappa_exact[n], " num factor:", sp.factor(kappa_exact[n]))
