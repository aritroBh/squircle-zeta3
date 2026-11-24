# s03: exact Beta-moment table J(n,k,a,b) = int_0^1 t^n (1-t)^k log^a t log^b(1-t) dt, a+b <= 3,
# as rational + rational*z2 + rational*z3 (z2 = zeta(2), z3 = zeta(3)).  Euler gamma cancels identically.
# Proof: log B(x+u, y+v) = sum_j [psi^{(j-1)}(x)u^j + psi^{(j-1)}(y)v^j - psi^{(j-1)}(x+y)(u+v)^j]/j!,
# psi(n) = -g + H_{n-1}, psi'(n) = z2 - H^{(2)}_{n-1}, psi''(n) = -2 z3 + 2 H^{(3)}_{n-1}.
import sympy as sp
from functools import lru_cache
u, v, z2, z3, g = sp.symbols('u v z2 z3 g')
def H(n, r): return sum(sp.Rational(1, j**r) for j in range(1, n))
def lg_series(n):  # lnGamma(n + w) - lnGamma(n) as polynomial coefficients [w^1, w^2, w^3]
    return [-g + H(n, 1), (z2 - H(n, 2))/2, (-2*z3 + 2*H(n, 3))/6]
@lru_cache(maxsize=None)
def J_table(n, k):
    x, y = n + 1, k + 1
    a, b, c = lg_series(x), lg_series(y), lg_series(x + y)
    S = sum(a[j]*u**(j+1) + b[j]*v**(j+1) - c[j]*(u+v)**(j+1) for j in range(3))
    S = sp.expand(S)
    E = sp.Integer(1); P = sp.Integer(1)
    for j in range(1, 4):
        P = sp.expand(P*S)
        E += P/sp.factorial(j)
    E = sp.expand(E)
    B0 = sp.factorial(n)*sp.factorial(k)/sp.factorial(n + k + 1)
    tab = {}
    for A in range(4):
        for Bb in range(4 - A):
            co = E.coeff(u, A).coeff(v, Bb)
            co = sp.expand(co.subs({u: 0, v: 0}))
            tab[(A, Bb)] = sp.expand(B0*sp.factorial(A)*sp.factorial(Bb)*co)
    return tab
def J(n, k, a, b): return J_table(n, k)[(a, b)]
if __name__ == '__main__':
    import mpmath as mp
    mp.mp.dps = 40
    worst = 0
    for n in range(0, 5):
        for k in range(0, 3):
            for a in range(4):
                for b in range(4 - a):
                    ex = J(n, k, a, b)
                    assert g not in ex.free_symbols, (n, k, a, b, ex)
                    val = ex.subs({z2: sp.pi**2/6, z3: sp.zeta(3)})
                    num = mp.quad(lambda t: t**n*(1-t)**k*mp.log(t)**a*mp.log(1-t)**b, [0, 0.5, 1])
                    worst = max(worst, abs(sp.N(val, 40) - num))
                    z3c = ex.coeff(z3)
                    if (a, b) in [(2, 1), (1, 2)] and k == 0:
                        assert z3c == 2*sp.Rational(1, n+1), (n, a, b, z3c)
                    elif (a, b) not in [(2, 1), (1, 2)]:
                        assert z3c == 0, (n, k, a, b, z3c)
    print("Beta table: Euler gamma cancels; numeric max error", worst)
    print("z3-coefficient rule verified: 2/(n+1) for (2,1),(1,2) [k=0], 0 for all other (a,b), a+b<=3")
    print("examples: J(0,0,2,1) =", J(0, 0, 2, 1), "; J(0,0,3,0) =", J(0, 0, 3, 0), "; J(2,1,1,1) =", J(2, 1, 1, 1))
