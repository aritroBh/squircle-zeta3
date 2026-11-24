"""
F1 step 4: the circle end p = 2 through the support function and Parseval.

eps = 2/p - 1 (as in research/05_final_verify/circle_end), so p = 2/(1+eps) and the dual exponent is q = 2/(1-eps):
polar duality is eps -> -eps.  With C = cos^2, S = sin^2, a = log|cos|, b = log|sin| and the cumulants of the
two-point law X = a (prob C), b (prob S):
    k1 = C a + S b,   k2 = C S (a-b)^2,   k3 = C S (S-C) (a-b)^3,
one gets EXACTLY (formal series)  log h_p = k1 eps + k2 eps^2 + (k2 + 2 k3/3) eps^3 + O(eps^4),
so h_p = 1 + k1 eps + (k2 + k1^2/2) eps^2 + (k2 + 2k3/3 + k1 k2 + k1^3/6) eps^3 + ...
(1) P = int_0^{2pi} h_p : integrate every monomial exactly with Beta derivatives; compare with circle_end's P_1..P_3.
(2) Hurwitz/Parseval: P^2 - 4 pi A = 2 pi^2 sum_{k>=1} (16k^2-1) a_k(h_p)^2 ; expand to eps^3:
       a_k = alpha_k eps + beta_k eps^2 + ..., alpha_k = [k1]_k = 1/(2k(4k^2-1)),
       beta_k = [k2 + k1^2/2]_k  (explicit, contains harmonic numbers H_k and odd harmonic numbers O_{2k}).
    eps^2:  2 pi^2 sum (16k^2-1) alpha_k^2      = D_2   (only pi^2 and rationals: the parity lock)
    eps^3:  4 pi^2 sum (16k^2-1) alpha_k beta_k = D_3   (contains zeta(3))
"""
import sympy as sp
import mpmath as mp
from math import factorial

# ---------- (1) exact perimeter coefficients from the support function ----------
eps = sp.Symbol('eps')
C, a, b = sp.symbols('C a b')
S = 1 - C
k1 = C * a + S * b
k2 = C * S * (a - b) ** 2
k3 = C * S * (S - C) * (a - b) ** 3

# check the cumulant expansion of log h symbolically against the definition, numerically at a point
def logh_exact(e, th):
    q = 2 / (1 - e)
    return mp.log(abs(mp.cos(th)) ** q + abs(mp.sin(th)) ** q) / q
mp.mp.dps = 40
th0 = mp.mpf('0.41')
vals = {C: mp.cos(th0) ** 2, a: mp.log(mp.cos(th0)), b: mp.log(mp.sin(th0))}
ser = [sp.N(k1.subs(vals), 40), sp.N(k2.subs(vals), 40), sp.N((k2 + sp.Rational(2, 3) * k3).subs(vals), 40)]
taylor = mp.taylor(lambda e: logh_exact(e, th0), 0, 3)
print("log h Taylor coeffs (numeric):", [mp.nstr(t, 20) for t in taylor[1:]])
print("cumulant formula            :", [mp.nstr(mp.mpf(str(x)), 20) for x in ser])

h_series = sp.expand(1 + k1 * eps + (k2 + k1 ** 2 / 2) * eps ** 2 + (k2 + sp.Rational(2, 3) * k3 + k1 * k2 + k1 ** 3 / 6) * eps ** 3)

# Beta-derivative integration: J(i,j,m,n) = int_0^{pi/2} cos^{2i} sin^{2j} a^m b^n
s_, t_ = sp.symbols('s t')
def Jmono(i, j, m, n):
    x0, y0 = sp.Rational(2 * i + 1, 2), sp.Rational(2 * j + 1, 2)
    order = m + n
    expo = 0
    for kk in range(1, order + 1):
        expo += (sp.polygamma(kk - 1, x0) * s_ ** kk + sp.polygamma(kk - 1, y0) * t_ ** kk
                 - sp.polygamma(kk - 1, x0 + y0) * (s_ + t_) ** kk) / sp.factorial(kk)
    Bser = sp.series(sp.exp(expo), s_, 0, m + 1).removeO()
    Bser = sp.series(Bser, t_, 0, n + 1).removeO()
    coef = sp.expand(Bser).coeff(s_, m).coeff(t_, n) * sp.factorial(m) * sp.factorial(n)
    B0 = sp.gamma(x0) * sp.gamma(y0) / sp.gamma(x0 + y0)
    return sp.Rational(1, 2 ** (m + n)) * B0 * coef / 2

def integrate_poly(expr):
    expr = sp.expand(expr)
    poly = sp.Poly(expr, C, a, b)
    total = 0
    for (i, m, n), cf in poly.terms():
        # C^i (S expanded already since S = 1-C)
        total += cf * Jmono(i, 0, m, n)
    return sp.simplify(sp.expand_func(sp.expand(total)))

def pg_exact(kk, x):
    """exact polygamma at positive integers and half-integers"""
    x = sp.Rational(x)
    if x.q == 1:
        n = int(x)
        if kk == 0:
            return -sp.EulerGamma + sum(sp.Rational(1, j) for j in range(1, n))
        return (-1) ** (kk + 1) * sp.factorial(kk) * (sp.zeta(kk + 1) - sum(sp.Rational(1, j ** (kk + 1)) for j in range(1, n)))
    assert x.q == 2
    n = int(x - sp.Rational(1, 2))
    if kk == 0:
        return -sp.EulerGamma - 2 * sp.log(2) + sum(sp.Rational(2, 2 * j - 1) for j in range(1, n + 1))
    return (-1) ** (kk + 1) * sp.factorial(kk) * ((2 ** (kk + 1) - 1) * sp.zeta(kk + 1)
            - sum(sp.Rational(2, 2 * j - 1) ** (kk + 1) for j in range(1, n + 1)))
def red(e):
    e = e.replace(lambda w: isinstance(w, sp.polygamma), lambda w: pg_exact(int(w.args[0]), w.args[1]))
    e = sp.expand(sp.expand_log(sp.expand(e), force=True))
    return sp.nsimplify(e) if False else e
lam, z3 = sp.log(2), sp.zeta(3)
P_support = {}
for kk in range(0, 4):
    term = h_series.coeff(eps, kk)
    P_support[kk] = red(4 * integrate_poly(term) / sp.pi)
    print(f"P_{kk}/pi (support function) =", P_support[kk])
l2, Z3, PI = sp.log(2), sp.zeta(3), sp.pi
P_circle_end = {0: sp.Integer(2), 1: 1 - 2 * l2, 2: (96 * l2 ** 2 - 96 * l2 + 11 * PI ** 2 - 72) / 96,
                3: (-128 * l2 ** 3 + 192 * l2 ** 2 - 44 * l2 * PI ** 2 + 288 * l2 + 37 * PI ** 2 - 276 * Z3 + 24) / 384}
for kk in range(4):
    d = sp.simplify(red(P_support[kk] - P_circle_end[kk]))
    print(f"P_{kk}: support - circle_end(arclength) = {d}   (numeric {sp.N(P_support[kk]-P_circle_end[kk], 30)})")

# ---------- (2) Hurwitz/Parseval deficit ----------
# log A in eps: A = 4 Gamma((3+eps)/2)^2/Gamma(2+eps)
logA = sp.log(4) + 2 * sp.loggamma((3 + eps) / 2) - sp.loggamma(2 + eps)
lA = [red(sp.diff(logA, eps, n).subs(eps, 0) / sp.factorial(n)) for n in range(4)]
print("\nlog A Taylor coefficients:", lA)
A_ser = sp.series(sp.exp(sum(lA[n] * eps ** n for n in range(4))), eps, 0, 4).removeO()
P_ser = sp.pi * sum(P_circle_end[n] * eps ** n for n in range(4))
D_ser = sp.expand(sp.series(P_ser ** 2 - 4 * sp.pi * A_ser, eps, 0, 4).removeO())
D = [red(D_ser.coeff(eps, n)) for n in range(4)]
print("Deficit D = P^2 - 4 pi A: D_0, D_1 =", D[0], ",", D[1])
print("D_2 =", D[2], "   D_2/(2 pi^2) =", red(D[2] / (2 * sp.pi ** 2)))
print("D_3 =", sp.simplify(D[3]))
D3n = red(D[3] / (4 * sp.pi ** 2))
print("D_3/(4 pi^2) =", D3n, " =", sp.N(D3n, 40))

mp.mp.dps = 50
lamn = mp.log(2)
# epsilon^2 Parseval sum (rational summand)
k = sp.Symbol('k', integer=True, positive=True)
summand2 = (16 * k ** 2 - 1) / (4 * k ** 2 * (4 * k ** 2 - 1) ** 2)
S2 = sp.summation(summand2, (k, 1, sp.oo))
print("\nsum (16k^2-1) alpha_k^2 =", sp.simplify(S2), "  vs D_2/(2pi^2) =", red(D[2] / (2 * sp.pi ** 2)),
      " diff", sp.simplify(red(S2 - D[2] / (2 * sp.pi ** 2))))

# explicit beta_k
def H(n):
    return mp.fsum(mp.mpf(1) / i for i in range(1, n + 1))
def beta_list(K):
    """beta_k for k = 1..K via explicit Fourier convolution formulas (O(K))"""
    Hk = [mp.mpf(0)] * (K + 3); Ok = [mp.mpf(0)] * (K + 3)   # Ok[j] = O_{2j} = sum_{i=1}^j 1/(2i-1)
    for j in range(1, K + 3):
        Hk[j] = Hk[j - 1] + mp.mpf(1) / j
        Ok[j] = Ok[j - 1] + mp.mpf(1) / (2 * j - 1)
    Y2 = lambda j: (mp.pi ** 2 / 2) if j == 0 else 4 * Ok[j] / j          # a_j of Y^2 (a_0 convention)
    g = lambda m: mp.mpf(-2) if m == 0 else -mp.mpf(4 * m) / (4 * m * m - 1)  # a_m of cos(2th) Y
    psi = mp.digamma; gam = mp.euler
    Sfun = lambda c: (psi(1 + mp.mpf(c) / 2) + gam) / c                     # sum_{n>=1} 1/(n(2n+c))
    Tfun = lambda kk, e: (psi(1 + kk) - psi(1 + mp.mpf(e) / 2)) / (2 * kk - e)  # sum_{m>=1} 1/((m+k)(2m+e))
    out = []
    for kk in range(1, K + 1):
        sum_up = -(Sfun(2 * kk - 1) + Sfun(2 * kk + 1))                           # sum_n (1/n) g_{n+k}
        fin = -((Hk[kk - 1] + 2 * Ok[kk - 1]) / (2 * kk - 1) + (Hk[kk - 1] + 2 * (Ok[kk] - 1)) / (2 * kk + 1))  # n=1..k-1
        mid = g(0) / kk
        tail = -(Tfun(kk, -1) + Tfun(kk, 1))                                       # n>k
        XcY = -(sum_up + fin + mid + tail) / 2
        X2 = (2 * Hk[kk] - mp.mpf(1) / kk) / kk
        beta = (X2 / 8 + lamn / (2 * kk) - XcY / 4 + lamn / 2 * g(kk)
                + mp.mpf(3) / 16 * Y2(kk) - (Y2(kk - 1) + Y2(kk + 1)) / 32)
        out.append(beta)
    return out

# check beta_k against direct quadrature of k2 + k1^2/2
def F2(th):
    Cn = mp.cos(th) ** 2; an = mp.log(abs(mp.cos(th))); bn = mp.log(abs(mp.sin(th)))
    K1 = Cn * an + (1 - Cn) * bn; K2 = Cn * (1 - Cn) * (an - bn) ** 2
    return K2 + K1 ** 2 / 2
def coeff(f, kk):
    return 4 / mp.pi * mp.quad(lambda t: f(t) * mp.cos(4 * kk * t), [0, mp.pi / 8, mp.pi / 4, 3 * mp.pi / 8, mp.pi / 2])
bl = beta_list(6)
for kk in range(1, 7):
    print(f"beta_{kk}: formula {mp.nstr(bl[kk-1], 25)}  quad diff {mp.nstr(bl[kk-1] - coeff(F2, kk), 3)}")
# also alpha_k
K1f = lambda th: mp.cos(th) ** 2 * mp.log(abs(mp.cos(th))) + mp.sin(th) ** 2 * mp.log(abs(mp.sin(th)))
for kk in [1, 2, 7]:
    print(f"alpha_{kk}: quad {mp.nstr(coeff(K1f, kk), 20)} vs 1/(2k(4k^2-1)) = {mp.nstr(mp.mpf(1)/(2*kk*(4*kk*kk-1)), 20)}")

# the eps^3 Parseval sum with tail extrapolation
Kmax = 16000
bl = beta_list(Kmax)
terms = [(16 * kk * kk - 1) * bl[kk - 1] / (2 * kk * (4 * kk * kk - 1)) for kk in range(1, Kmax + 1)]
partial = []
acc = mp.mpf(0)
for kk, tt in enumerate(terms, 1):
    acc += tt
    partial.append(acc)
# fit S(K) = S + (A log K + B)/K^3 + (C log K + D)/K^4 + (E log K + F)/K^5 using several K
Ks = [4000, 6000, 8000, 10000, 12000, 14000, 16000]
def basis_row(KK):
    lk = mp.log(KK)
    return [1, lk / KK ** 3, mp.mpf(1) / KK ** 3, lk / KK ** 4, mp.mpf(1) / KK ** 4, lk / KK ** 5, mp.mpf(1) / KK ** 5]
scale = [abs(x) for x in basis_row(Ks[0])]
M = mp.matrix(len(Ks), 7); rhs = mp.matrix(len(Ks), 1)
for r, KK in enumerate(Ks):
    row = basis_row(KK)
    for cidx in range(7):
        M[r, cidx] = row[cidx] / scale[cidx]
    rhs[r] = partial[KK - 1]
sol = mp.lu_solve(M, rhs)
# a second, smaller fit as a stability check
Ks2 = [6000, 9000, 12000, 16000, 14000]
M2 = mp.matrix(5, 5); r2 = mp.matrix(5, 1)
for r, KK in enumerate(Ks2):
    row = basis_row(KK)[:5]
    for cidx in range(5):
        M2[r, cidx] = row[cidx] / scale[cidx]
    r2[r] = partial[KK - 1]
sol2 = mp.lu_solve(M2, r2)
print("tail-fit stability: 7-term vs 5-term extrapolation differ by", mp.nstr(sol[0] - sol2[0], 3))
S3 = sol[0]
target = mp.mpf(str(sp.N(D3n, 45)))
print(f"\nsum_(k<=16000) (16k^2-1) alpha_k beta_k = {mp.nstr(partial[-1], 25)}")
print(f"extrapolated sum                        = {mp.nstr(S3, 25)}")
print(f"D_3/(4 pi^2) exact                      = {mp.nstr(target, 25)}   diff {mp.nstr(S3 - target, 3)}")
zc = sp.expand(D3n).coeff(sp.zeta(3))
print("coefficient of zeta(3) in D_3/(4pi^2):", zc)
