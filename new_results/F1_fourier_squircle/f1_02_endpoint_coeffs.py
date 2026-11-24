"""
F1 step 2: Fourier coefficients of the squircle functions at the endpoints p = 1, 2, infinity,
and their first p-derivatives (question (a)).

Convention: every function below has the symmetry of the square, so it is even, pi/2-periodic and
symmetric about pi/4. We write f = a_0/2 + sum_{k>=1} a_k cos(4k theta), a_k = (4/pi) int_0^{pi/2} f cos(4k theta).

Closed forms claimed (proofs in result.md):
 (S1) support function h_1 = max(|c|,|s|) (diamond):   a_k = -(-1)^k 4 sqrt2 / (pi (16k^2-1))
 (S2) support function h_inf = |c|+|s| (square):        a_k = -8 / (pi (16k^2-1))
 (R1) radial function rho_1 = 1/(|c|+|s|):             a_k = (8/pi)(-1)^k T_{4k},  T_N = sum_{n odd > N} chi_8(n)/n,
                                                         a_0 = 4 sqrt2 L / pi   (L = log(1+sqrt2))
 (G1) log rho_1 = -log(|c|+|s|):                       a_k = (-1)^k (2/(pi k)) sum_{j>=k} (-1)^j/(2j+1)  (k>=1),
                                                         a_0 = -log 2 + 4G/pi... (checked numerically below)
 (U)  u = d/dp log rho_p |_{p=2} = d/dp h_p |_{p=2} = -(c^2 log|c| + s^2 log|s|)/2:
                                                         a_k = -1/(4k(4k^2-1)),  a_0 = log 2 - 1/2
 (Gam) |c|^p+|s|^p = rho_p^{-p} = sum_m gamma_m(p) e^{4 i m theta},
          gamma_m(p) = 2^{1-p} Gamma(p+1) / (Gamma(1+p/2+2m) Gamma(1+p/2-2m))
 (V)  v = d/dp log rho_p |_{p=1}:  in u = tan theta, v = log(1+u) - u log u/(1+u); its coefficients are
      identified by PSLQ below on the basis {1, log2, pi, pi^2, pi log2, log^2 2, G} (after multiplying by pi/4).
"""
import mpmath as mp

mp.mp.dps = 70
pi = mp.pi; L = mp.log(1 + mp.sqrt(2)); lam = mp.log(2); G = mp.catalan; s2 = mp.sqrt(2)

def coeff(f, k, extra=None):
    """a_k = (4/pi) int_0^{pi/2} f cos(4k th) dth, split at pi/8, pi/4, 3pi/8"""
    pts = [0, pi / 8, pi / 4, 3 * pi / 8, pi / 2]
    return 4 / pi * mp.quad(lambda t: f(t) * mp.cos(4 * k * t), pts)

c = mp.cos; s = mp.sin
h1 = lambda t: max(abs(c(t)), abs(s(t)))
hinf = lambda t: abs(c(t)) + abs(s(t))
rho1 = lambda t: 1 / (abs(c(t)) + abs(s(t)))
logrho1 = lambda t: -mp.log(abs(c(t)) + abs(s(t)))
def uu(t):
    C, S = c(t), s(t)
    return -(C * C * mp.log(abs(C)) + S * S * mp.log(abs(S))) / 2

def chi8(n):
    return {1: 1, 3: -1, 5: -1, 7: 1}.get(n % 8, 0)

def T_tail(N):
    # sum_{n odd > N} chi8(n)/n = L/sqrt2 - sum_{n odd <= N} chi8(n)/n
    return L / s2 - mp.fsum(mp.mpf(chi8(n)) / n for n in range(1, N + 1, 2))

def leib_tail(k):
    return pi / 4 - mp.fsum(mp.mpf(-1) ** j / (2 * j + 1) for j in range(k))

print("== (S1),(S2),(R1),(G1),(U): numerical vs closed form, k = 0..6 ==")
worst = 0
for k in range(0, 7):
    rows = []
    a = coeff(h1, k); cf = -(-1) ** k * 4 * s2 / (pi * (16 * k * k - 1)); rows.append(('h1', a, cf))
    a = coeff(hinf, k); cf = -8 / (pi * (16 * k * k - 1)); rows.append(('hinf', a, cf))
    a = coeff(rho1, k); cf = (4 * s2 * L / pi) if k == 0 else 8 / pi * (-1) ** k * T_tail(4 * k); rows.append(('rho1', a, cf))
    if k >= 1:
        a = coeff(logrho1, k); cf = (-1) ** k * 2 / (pi * k) * leib_tail(k); rows.append(('logrho1', a, cf))
    a = coeff(uu, k); cf = (lam - mp.mpf(1) / 2) if k == 0 else -mp.mpf(1) / (4 * k * (4 * k * k - 1)); rows.append(('u', a, cf))
    for name, a, cf in rows:
        worst = max(worst, abs(a - cf))
        print(f"k={k:2d} {name:8s} {mp.nstr(a, 25):>32s}  diff {mp.nstr(a - cf, 3)}")
print("worst diff:", mp.nstr(worst, 3))

a0_logrho1 = coeff(logrho1, 0)
print("\nlog rho_1: a_0 =", mp.nstr(a0_logrho1, 40))
rel = mp.pslq([a0_logrho1, 1, lam, G / pi, L], maxcoeff=10 ** 6, maxsteps=10 ** 6)
print("  PSLQ [a0, 1, log2, G/pi, L]:", rel)
print("  check -log2 + 4G/pi... value of -lam + 4*G/pi - ? :", mp.nstr(a0_logrho1 - (-lam + 4 * G / pi), 20))

print("\n== (Gam): Fourier coefficients of |c|^p+|s|^p  (complex form, frequency 4m) ==")
def gam(m, p):
    return 2 ** (1 - p) * mp.gamma(p + 1) * mp.rgamma(1 + p / 2 + 2 * m) * mp.rgamma(1 + p / 2 - 2 * m)
for p in [mp.mpf('0.7'), mp.mpf('1.3'), mp.mpf(3), mp.mpf('5.5')]:
    for m in [0, 1, 2, 3]:
        num = coeff(lambda t: abs(c(t)) ** p + abs(s(t)) ** p, m) / (1 if m else 1)
        # complex coefficient = a_m/2 for m>=1, a_0/2 for m=0
        num = num / 2
        print(f"p={mp.nstr(p,3):4s} m={m}  quad {mp.nstr(num, 22):>28s}  gamma-formula diff {mp.nstr(num - gam(m, p), 3)}")

print("\n== p-derivatives of gamma_m at p=2 (rational) and p=1 (log 2 + odd harmonic sums) ==")
for m in range(0, 5):
    d2 = mp.diff(lambda p: gam(m, p), 2)
    cf2 = (mp.mpf(1) / 2 - lam) if m == 0 else mp.mpf(1) / (4 * m * (4 * m * m - 1))
    d1 = mp.diff(lambda p: gam(m, p), 1)
    Om = lambda N: mp.fsum(mp.mpf(1) / (2 * j - 1) for j in range(1, N + 1))
    if m == 0:
        cf1 = 4 / pi * (lam - 1)
    else:
        cf1 = -4 / (pi * (16 * m * m - 1)) * (lam + 1 - 2 * Om(2 * m - 1) - mp.mpf(8 * m) / (16 * m * m - 1))
    print(f"m={m}: d/dp gamma_m(2) = {mp.nstr(d2, 25):>30s}  vs {mp.nstr(cf2,12)} diff {mp.nstr(d2-cf2,3)}",
          f"\n      d/dp gamma_m(1) = {mp.nstr(d1, 25):>30s}  vs formula diff {mp.nstr(d1 - cf1, 3)}")

print("\n== (V): coefficients of v = d/dp log rho_p at p = 1 ==")
def v_theta(t):
    C, S = abs(c(t)), abs(s(t))
    return mp.log(C + S) - (C * mp.log(C) + S * mp.log(S)) / (C + S)
# cross-check v by finite differences in p
def logrho(p, t):
    return -mp.log(abs(c(t)) ** p + abs(s(t)) ** p) / p
t0 = mp.mpf('0.37')
print("v(t0) direct", mp.nstr(v_theta(t0), 25), " finite-diff", mp.nstr(mp.diff(lambda p: logrho(p, t0), 1), 25))
basis_names = ['1', 'log2', 'pi', 'pi^2', 'pi*log2', 'log2^2', 'G']
basis = [1, lam, pi, pi ** 2, pi * lam, lam ** 2, G]
for k in range(0, 5):
    vk = coeff(v_theta, k)
    x = vk * pi / 4
    rel = mp.pslq([x] + basis, maxcoeff=10 ** 8, maxsteps=10 ** 7)
    if rel:
        res = abs(x * rel[0] + mp.fsum(r * b for r, b in zip(rel[1:], basis)))
        terms = " + ".join(f"({mp.mpf(-r)/rel[0]})*{n}" for r, n in zip(rel[1:], basis_names) if r)
        print(f"k={k}: a_k(v) = {mp.nstr(vk, 30)};  (pi/4) a_k = {terms};  residual {mp.nstr(res, 3)}, max|coef| {max(abs(r) for r in rel)}")
    else:
        print(f"k={k}: a_k(v) = {mp.nstr(vk, 30)}; no relation")
