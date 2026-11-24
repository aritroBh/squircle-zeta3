import sympy as sp
from sympy import symbols, Rational, factorial, binomial, sqrt, exp, series, simplify

NMAX = 12
Z = {n: symbols(f'Z{n}') for n in range(2, NMAX+1)}
eta, mu, s, d, t, q = symbols('eta mu s d t q')

def ts_mul(A, B, N=NMAX):
    R = {}
    for i, aa in A.items():
        for j, bb in B.items():
            if i + j <= N:
                R[i+j] = R.get(i+j, 0) + aa*bb
    return {k: sp.expand(v) for k, v in R.items() if v != 0}

def ts_exp(A, N=NMAX):
    R = {0: sp.Integer(1)}
    for n in range(1, N+1):
        acc = 0
        for k in range(1, n+1):
            if k in A and (n-k) in R:
                acc += k * A[k] * R[n-k]
        val = sp.expand(acc / n)
        if val != 0: R[n] = val
    return R

# even part of sqrt(cosh z + q sinh z): c_k(q) polynomial in q
# g(z)^2 = cosh z + q sinh z; use g_e via series in z of sqrt, then pick even part coefficients
z = symbols('z')
gz = sqrt(sp.cosh(z) + q*sp.sinh(z))
gser = series(gz, z, 0, NMAX+1).removeO()
ge = sum((gser.coeff(z, k) + (-1)**k * gser.coeff(z, k))/2 * z**k for k in range(0, NMAX+1))
cq = {k: sp.expand(ge.coeff(z, 2*k)) for k in range(NMAX//2+1)}
print("c_k(q):", {k: cq[k] for k in range(4)})

mA = {n: (-1)**(n+1) * Z[n] / n for n in range(2, NMAX+1)}
K = ts_exp(mA)

def bs_exp(A, N=NMAX):
    R = {(0,0): sp.Integer(1)}
    for tot in range(1, N+1):
        for i in range(0, tot+1):
            j = tot - i
            if i > 0:
                acc = 0
                for (i1,j1), aa in A.items():
                    if i1 >= 1 and (i-i1, j-j1) in R:
                        acc += i1 * aa * R[(i-i1, j-j1)]
                val = sp.expand(acc / i)
            else:
                acc = 0
                for (i1,j1), aa in A.items():
                    if i1 == 0 and j1 >= 1 and (0, j-j1) in R:
                        acc += j1 * aa * R[(0, j-j1)]
                val = sp.expand(acc / j) if j > 0 else sp.Integer(1)
            if val != 0: R[(i,j)] = val
    return R

Bser = {}
for n in range(2, NMAX+1):
    for r in range(0, n+1, 2):
        Bser[(n-r, r)] = Bser.get((n-r, r), 0) + (-1)**n * Z[n] * 2 * binomial(n, r) / n
Bser = {k: sp.expand(v) for k, v in Bser.items() if v != 0}
Hser = bs_exp(Bser)
Hder = {}
for k in range(NMAX//2+1):
    Hd = {}
    for i in range(0, NMAX+1):
        v = Hser.get((i, 2*k), 0)
        if v != 0: Hd[i] = sp.expand(factorial(2*k) * v)
    Hder[k] = Hd

Sq = {}
for k in range(NMAX//2+1):
    for i, v in Hder[k].items():
        Sq[2*k+i] = Sq.get(2*k+i, 0) + cq[k] * v / 2**i
Sq = ts_mul(K, Sq)

# kappa_n(q) exact from S_q: kappa = 4(-1)^n [Z_n eta^n] S_q
A_, B_ = (1+q)/2, (1-q)/2
Rqt = (sqrt(A_*exp(4*t) + B_) + sqrt(A_ + B_*exp(4*t)))/2
Rqser = series(Rqt, t, 0, NMAX+1).removeO()
print(f"\n{'n':>3} {'kappa_n(q) from S_q':>40} {'closed form':>40} match")
for n in range(3, 12, 2):
    kq_exact = sp.expand(Sq.get(n, 0)).coeff(Z[n], 1) * 4 * (-1)**n
    Vq = sp.expand(factorial(n) * Rqser.coeff(t, n))
    kq_form = sp.expand(Rational(4*(-1)**n, n) * (1 - Rational(2)**(1-n) * Vq))
    print(f"{n:>3} {str(kq_exact):>40} {str(kq_form):>40} {sp.simplify(kq_exact-kq_form)==0}")

# in terms of a,b: A = a^2/(a^2+b^2), B = b^2/(a^2+b^2); V_3(q) explicit
a2, b2 = symbols('a2 b2', positive=True)
V3 = sp.expand(factorial(3)*Rqser.coeff(t,3))
print("\nV_3(q) =", sp.factor(V3), "  kappa_3(q) =", sp.factor(sp.expand(Rational(-4,3)*(1-Rational(1,4)*V3))))
