"""extra_checks.py -- checks behind the statements added or changed in revision round 1.

Independent of the other scripts (imports nothing from the other folders).
Run:  python extra_checks.py
Sections:
  1. Faulhaber's formula with (-1)^j B_j and B_1 = -1/2 (Section V-A), exact.
  2. omega_n = (4/n)(2^(1-n) V_n - 1), the coefficient of zeta(n) in a_n (Prop. a_n), odd n <= 21,
     against an independent symbolic expansion with formal zeta symbols for n <= 9.
  3. Ahmed's Table IV at n = 10^6 from the exact formula of Theorem P2; the gap C(theta) - C_n-hat and
     the second-order term that explains it.
  4. The remainder O(log n / n^3) in Theorem P2asym (exact formula vs the truncated expansion).
  5. The smoothed corner of interior angle alpha: length of the oblique quarter-squircle = 2 + 2u M(u; alpha).
  6. Li_3(1/2) = 7/8 zeta(3) - pi^2/12 log 2 + log^3 2 / 6 (used in Section IX-C).
"""
import sys
import time
import sympy as sp
from mpmath import mp, mpf, quad, hyp2f1, legenp, sin, cos, log, pi, sqrt, polylog, zeta, nstr, li, inf

LOG = []
T0 = time.time()


def check(tag, ok, detail=""):
    LOG.append((tag, ok))
    print(f"{'PASS' if ok else 'FAIL'} | {tag} {detail}", flush=True)


# ---------------------------------------------------------------- 1. Faulhaber
print("== 1. Faulhaber's formula (Section V-A) ==")
B = [sp.Rational(1)]
for n in range(1, 12):  # sum_{k=0}^{n} C(n+1,k) B_k = 0  gives B_1 = -1/2
    B.append(-sum(sp.binomial(n + 1, k) * B[k] for k in range(n)) / (n + 1))
assert B[1] == sp.Rational(-1, 2) and B[2] == sp.Rational(1, 6) and B[4] == sp.Rational(-1, 30)
ok_new, ok_old_fails = True, True
for m in range(1, 9):
    for N in range(1, 13):
        lhs = sum(sp.Integer(k) ** m for k in range(1, N + 1))
        new = sp.Rational(1, m + 1) * sum(sp.binomial(m + 1, j) * (-1) ** j * B[j] * sp.Integer(N) ** (m + 1 - j)
                                          for j in range(m + 1))
        old = sp.Rational(1, m + 1) * sum(sp.binomial(m + 1, j) * B[j] * sp.Integer(N) ** (m + 1 - j)
                                          for j in range(m + 1))
        ok_new &= (new == lhs)
        ok_old_fails &= (old == sum(sp.Integer(k) ** m for k in range(0, N)))
check("sum_{k=1}^N k^m = 1/(m+1) sum_j C(m+1,j) (-1)^j B_j N^(m+1-j), m<=8, N<=12 (exact)", ok_new)
check("the printed form with B_j (B_1=-1/2) equals sum_{k=0}^{N-1} k^m instead (exact)", ok_old_fails)

# ---------------------------------------------------------------- 2. omega_n
print("== 2. zeta(n) in a_n (Proposition a_n) ==")
t = sp.symbols("t")
NV = 21
ser = sp.series(sp.exp(t) * sp.sqrt(sp.cosh(2 * t)), t, 0, NV + 1).removeO()
V = {n: sp.factorial(n) * ser.coeff(t, n) for n in range(NV + 1)}
chk = (sp.exp(t) * sp.sqrt(sp.cosh(2 * t)) - sp.sqrt((1 + sp.exp(4 * t)) / 2)).subs(t, sp.Rational(1, 3)).evalf(40)
check("e^t sqrt(cosh 2t) = sqrt((1+e^{4t})/2), so V_n = (d/dt)^n sqrt((1+e^{4t})/2) at 0", abs(chk) < 1e-35)
omega = {n: sp.Rational(4, n) * (sp.Integer(2) ** (1 - n) * V[n] - 1) for n in range(3, NV + 1, 2)}
print("   V_n (odd n):", {n: V[n] for n in range(3, NV + 1, 2)})
print("   omega_n    :", omega)
check("omega_3, omega_5, omega_7 = 1, -3/4, 129/16 (the zeta(n) coefficients printed in a_3, a_5, eq. (a7))",
      (omega[3], omega[5], omega[7]) == (1, sp.Rational(-3, 4), sp.Rational(129, 16)))
check("omega_n != 0 for every odd n <= 21", all(omega[n] != 0 for n in omega))
# independent: formal expansion with zeta symbols z_m, read off the coefficient of z_n in [d^n]
NS = 9
eta, s, d = sp.symbols("eta s d")
z = {m: sp.Symbol(f"z{m}") for m in range(2, NS + 1)}
eps = sp.symbols("eps")  # bookkeeping: total degree in (eta, s)
logA = sum((-1) ** m * z[m] / m * ((eps * eta / 2 + eps * s) ** m + (eps * eta / 2 - eps * s) ** m - (eps * eta) ** m)
           for m in range(2, NS + 1))
# only the part linear in the z's can carry z_n at degree n (products have degree >= n+2 when they contain z_n)
A_lin = 1 + logA
cosh_ser = sp.series(sp.sqrt(sp.cosh(sp.Symbol("x"))), sp.Symbol("x"), 0, NS + 2).removeO()
cks = {k: cosh_ser.coeff(sp.Symbol("x"), 2 * k) for k in range(0, NS // 2 + 1)}
expr = 0
for k in range(0, NS // 2 + 1):
    der = sp.diff(A_lin, s, 2 * k).subs(s, 0)
    expr += cks[k] * (eps * eta) ** (2 * k) * der
expr = sp.expand(expr)
ok = True
for n in range(3, NS + 1, 2):
    cz = sp.expand(expr).coeff(z[n])                    # polynomial in eps*eta
    cz = sp.expand(cz.subs(eps, 1))
    # P/sqrt2 = 4 * (P/(4 sqrt2)); eta = -d/(1+d); coefficient of d^n
    poly = sp.series(4 * cz.subs(eta, -d / (1 + d)), d, 0, n + 1).removeO()
    coeff_dn = sp.nsimplify(poly.coeff(d, n))
    lower = [sp.simplify(poly.coeff(d, j)) for j in range(0, n)]
    ok &= (coeff_dn == omega[n]) and all(c == 0 for c in lower)
    print(f"   n = {n}: coefficient of zeta({n}) in a_{n} (formal expansion) = {coeff_dn}; "
          f"zeta({n}) absent from a_0..a_{n-1}: {all(c == 0 for c in lower)}")
check("formal expansion reproduces omega_n for odd n <= 9 and zeta(n) is absent below d^n", ok)

# ---------------------------------------------------------------- 3. Ahmed's Table IV
print("== 3. Ahmed's Table IV at n = 10^6 (Theorem P2) ==")
mp.dps = 40


def Lam_eps(e, th):
    zz = sin(th / 2) ** 2
    return 2 * zz * hyp2f1(-e, 1 + e, 2, zz)


def Lam_m(m, th):
    c = cos(th)
    return (legenp(m - 1, 0, c) - legenp(m + 1, 0, c)) / (2 * m + 1)


def L_exact(n, th, M=60):
    e = 1 / mpf(n - 1)
    val = 1 + cos(th) + (pi * e / sin(pi * e)) * Lam_eps(e, th) * mpf(n) ** (-e)
    for m in range(1, M):
        val += (-1) ** m * e / (m - e) * Lam_m(m, th) * mpf(n) ** (-m)
    return val


def a1(th):
    w = 1 / mp.tan(th / 2) ** 2
    return -1 - w * log(cos(th / 2) ** 2)


def a2(th):
    zz = sin(th / 2) ** 2
    w = 1 / mp.tan(th / 2) ** 2
    u = w * log(1 - zz)
    return zeta(2) + mpf(1) / 2 - polylog(2, zz) - u ** 2 / 2


def c2(th):
    zz = sin(th / 2) ** 2
    return 2 * zz * (a2(th) - 1 + zz)


ahmed = {"pi/6": ("0.0045910", "0.0045781"), "pi/4": ("0.0225777", "0.0225492"),
         "pi/3": ("0.0684769", "0.0684278"), "pi/2": ("0.3068528", "0.3067523"),
         "2pi/3": ("0.8068528", "0.8066971"), "3pi/4": ("1.1444313", "1.1442509"),
         "5pi/6": ("1.5038583", "1.5036576")}
angle = {"pi/6": pi / 6, "pi/4": pi / 4, "pi/3": pi / 3, "pi/2": pi / 2, "2pi/3": 2 * pi / 3,
         "3pi/4": 3 * pi / 4, "5pi/6": 5 * pi / 6}
n = 10 ** 6
all_ok, gaps, pred_err = True, [], []
for key, th in angle.items():
    zz = sin(th / 2) ** 2
    C = 2 * zz + 4 * cos(th / 2) ** 2 * log(cos(th / 2))
    Chat = (n - 1) * (2 - L_exact(n, th)) - 2 * zz * log(n)
    pred = C - (zz * (log(n) - a1(th)) ** 2 + c2(th)) / (n - 1)
    okC = abs(C - mpf(ahmed[key][0])) < 5e-8
    okH = abs(Chat - mpf(ahmed[key][1])) < 5e-8
    all_ok &= okC and okH
    gaps.append(C - Chat)
    pred_err.append(abs(pred - Chat))
    print(f"   {key:>5}: C = {nstr(C, 10)} (Ahmed {ahmed[key][0]}), C_n-hat = {nstr(Chat, 12)} (Ahmed {ahmed[key][1]}), "
          f"gap = {nstr(C - Chat, 4)}, |2nd-order prediction - C_n-hat| = {nstr(abs(pred - Chat), 3)}")
check("C(theta) = -2 z a_1 and the exact C_n-hat match all 7 rows of Ahmed's Table IV to every printed digit", all_ok)
check(f"gap C - C_n-hat between {nstr(min(gaps), 3)} and {nstr(max(gaps), 3)}", 1.2e-5 < min(gaps) and max(gaps) < 2.1e-4)
check("second-order term -[z(log n - a_1)^2 + c_2(theta)]/(n-1) explains the gap to within 1e-9", max(pred_err) < 1e-9,
      f"(worst {nstr(max(pred_err), 3)})")

# ---------------------------------------------------------------- 4. remainder in Theorem P2asym
print("== 4. Remainder of Theorem P2asym ==")
mp.dps = 60
ok = True
for th in (pi / 3, pi / 2, 3 * pi / 4):
    zz = sin(th / 2) ** 2
    mu = mp.e ** a1(th)
    row = []
    for n in (10 ** 4, 10 ** 6, 10 ** 8):
        approx = 1 + cos(th) + 2 * zz * (mu / n) ** (1 / mpf(n - 1)) + c2(th) / mpf(n) ** 2
        r = (L_exact(n, th) - approx) * mpf(n) ** 3 / log(n)
        row.append(r)
    ok &= all(abs(x) < 10 for x in row)
    print(f"   theta = {nstr(th, 6)}: (L - truncation) n^3 / log n at n = 1e4, 1e6, 1e8: {[nstr(x, 6) for x in row]}")
check("(L_theta(n) - truncation after c_2/n^2) * n^3/log n stays bounded", ok)
mp.dps = 40
check("mu_{pi/2} = e^{a_1(pi/2)} = 2/e and c_2(pi/2) = pi^2/12",
      abs(mp.e ** a1(pi / 2) - 2 / mp.e) < mpf(10) ** -35 and abs(c2(pi / 2) - pi ** 2 / 12) < mpf(10) ** -35)

# ---------------------------------------------------------------- 5. smoothed corner
print("== 5. Smoothed corner of interior angle alpha (Section XII-B) ==")
mp.dps = 30


def corner_length(p, al):
    # gamma(t) = x e1 + y e2, x = t^(1/p), y = (1-t)^(1/p), e1.e2 = cos(al).  By symmetry, twice the
    # half t in [0, 1/2]; with t = w^p the endpoint singularity t^(u-1) disappears (u t^(u-1) dt = dw).
    u = 1 / mpf(p)
    def f(w):
        tt = w ** p
        r = (tt / (1 - tt)) ** (1 - u)          # = b/a
        return sqrt(1 + r ** 2 - 2 * cos(al) * r)
    return 2 * quad(f, [0, mpf(2) ** (-u)])


def M_corner(u, al):
    f = lambda g: ((1 + g) ** (-1 - u) * sqrt(1 + g ** (2 - 2 * u) - 2 * cos(al) * g ** (1 - u)) - 1) * g ** (u - 1)
    return quad(f, [0, 1])


ok = True
for p, al in ((5, 2 * pi / 3), (5, pi / 2), (12, pi / 6)):
    u = 1 / mpf(p)
    lhs, rhs = corner_length(p, al), 2 + 2 * u * M_corner(u, al)
    ok &= abs(lhs - rhs) < mpf(10) ** -25
    print(f"   p = {p}, alpha = {nstr(al, 6)}: length {nstr(lhs, 25)}  2+2uM = {nstr(rhs, 25)}")
# alpha = pi/2: quarter of the squircle perimeter, P/4 from eq. (1)
p = 5
# eq. (1) of the paper with the same substitution: P/4 = (1/p) int_0^1 sqrt(t^(2/p-2) + (1-t)^(2/p-2)) dt
P4 = 2 * quad(lambda w: sqrt(1 + ((w ** p) / (1 - w ** p)) ** (2 - 2 / mpf(p))), [0, mpf(2) ** (-1 / mpf(p))])
ok &= abs(P4 - corner_length(p, pi / 2)) < mpf(10) ** -25
print(f"   P(5)/4 from eq. (1) = {nstr(P4, 25)}")
check("length of x e1 + y e2 (x^p + y^p = 1, e1.e2 = cos alpha) = 2 + 2u M(u; alpha), u = 1/p; = P(p)/4 at alpha = pi/2", ok)

# ---------------------------------------------------------------- 6. Li3(1/2)
mp.dps = 50
v = polylog(3, mpf(1) / 2) - (mpf(7) / 8 * zeta(3) - pi ** 2 / 12 * log(2) + log(2) ** 3 / 6)
check("Li_3(1/2) = 7/8 zeta(3) - pi^2/12 log2 + log^3 2/6", abs(v) < mpf(10) ** -48, f"(resid {nstr(v, 3)})")

nf = sum(1 for _, ok in LOG if not ok)
print(f"\nSUMMARY: {len(LOG)} checks, {nf} failed, {time.time() - T0:.0f}s")
