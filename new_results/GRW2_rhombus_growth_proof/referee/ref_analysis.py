"""Referee checks built on In_direct.pkl (from ref_coeffs.py).

H_n: own series of h1 = sqrt(A^2+B^2 e^{-2z})/r, h2 = sqrt(B^2+A^2 e^{-2z})/r.
Checks: root test, E_n = I_n/n! - H_n decay (Lemma 4, with explicit C'_sigma), local constant c_1 at the
branch point by direct evaluation, Theorem B (four-point transfer, two terms) and its real form,
Lemma 6 mean square, Laplace representation (Lemma 1) at several eta.
"""
import os, pickle
from mpmath import (mp, mpf, mpc, sqrt, log, exp, pi, factorial, quad, inf, cos, sin, arg, fabs,
                    nstr, re, im)

mp.dps = 50
here = os.path.dirname(__file__)
data = pickle.load(open(os.path.join(here, "In_direct.pkl"), "rb"))
N = 40

def series_exp(c, N):  # coefficients of e^{c z}
    out, f = [], mpf(1)
    for k in range(N + 1):
        if k: f *= k
        out.append(c ** k / f)
    return out

def sqrt_series(c, N):
    s = [sqrt(c[0])]
    for n in range(1, N + 1):
        s.append((c[n] - sum(s[k] * s[n - k] for k in range(1, n))) / (2 * s[0]))
    return s

def Hcoeffs(A, B, N):
    r2 = A * A + B * B
    e = series_exp(mpf(-2), N)
    R1 = [B * B * x for x in e]; R1[0] += A * A
    R2 = [A * A * x for x in e]; R2[0] += B * B
    h1 = sqrt_series(R1, N); h2 = sqrt_series(R2, N)
    return [(x + y) / sqrt(r2) for x, y in zip(h1, h2)]

def Hfun(A, B, z):  # principal roots; valid in the strip |Im z| < pi/2 and on D(0,rho)
    r = sqrt(A * A + B * B)
    return (sqrt(A * A + B * B * exp(-2 * z)) + sqrt(B * B + A * A * exp(-2 * z))) / r

for (A, B) in [(1, 1), (1, 3), (2, 5)]:
    A, B = mpf(A), mpf(B)
    print("=" * 78)
    print(f"(A,B) = ({int(A)},{int(B)})")
    I = [mpf(x) for x in data[(int(A), int(B))]]
    r = sqrt(A * A + B * B); L = log(B / A); rho = sqrt(L * L + pi * pi / 4)
    th = arg(mpc(L, pi / 2))
    z = [mpc(L, pi / 2), mpc(L, -pi / 2), mpc(-L, pi / 2), mpc(-L, -pi / 2)]
    print(f"rho = {nstr(rho, 12)}, 1/rho = {nstr(1/rho, 12)}, theta0/pi = {nstr(th/pi, 8)}")
    print(f"I_0 = {nstr(I[0], 20)} (must be 1)")
    H = Hcoeffs(A, B, N)

    c = [(A / r) * sqrt(-2 * z[0]), (A / r) * sqrt(-2 * z[1]),
         (B / r) * sqrt(2 * z[0].conjugate()), (B / r) * sqrt(2 * z[0])]
    # local constants by direct evaluation along the radial segment t*z_j, t -> 1- (principal roots there)
    t = 1 - mpf(10) ** -14
    for j in range(4):
        zz = t * z[j]
        hj = (sqrt(A * A + B * B * exp(-2 * zz)) if j < 2 else sqrt(B * B + A * A * exp(-2 * zz))) / r
        print(f"c{j+1}: numeric {nstr(hj / sqrt(1 - t), 12)}  claimed {nstr(c[j], 12)}")
    if A == B:
        pass  # z3=z1, z4=z2: the formula below just adds the coincident contributions

    # sigma for Lemma 4 explicit check
    sig = mpf("0.9") * rho
    Cs = max(abs(Hfun(A, B, sig * exp(mpc(0, 2 * pi * k / 720)))) for k in range(720))
    ser = sum((2 * sig) ** j / factorial(j) ** 2 for j in range(1, 80))
    Cp = Cs * (1 + ser / 2)

    print(" n | |I_n/n!|^(1/n) | |H_n|^(1/n) | E_n(2rho)^n | Lemma4 ok | N_n | law1 | (I/n!-law2)rho^n n^3.5")
    rows = []
    for n in range(1, N + 1):
        x = I[n] / factorial(n)
        En = x - H[n]
        Nn = rho ** n * mpf(n) ** 1.5 * x
        beta = B * (-1) ** n * cos((n - mpf(1) / 2) * th) - A * sin((n - mpf(1) / 2) * th)
        law1 = -sqrt(2 * rho) / (r * sqrt(pi)) * beta
        four = sum(cj * zj ** (-n) * (1 + mpf(3) / (8 * n) - 3 * zj / (4 * n)) for cj, zj in zip(c, z))
        law2 = -(1 / (2 * sqrt(pi))) * mpf(n) ** -1.5 * four
        four1 = sum(cj * zj ** (-n) for cj, zj in zip(c, z))
        lead_real = re(-(1 / (2 * sqrt(pi))) * four1) * rho ** n
        assert fabs(lead_real - law1) < mpf(10) ** -30, (n, lead_real, law1)
        assert fabs(im(four)) < mpf(10) ** -30
        ok = fabs(En) <= Cp * (2 * sig) ** -n
        res = (x - re(law2)) * rho ** n * mpf(n) ** 3.5
        rows.append((n, Nn, law1))
        if n % 4 == 0 or n <= 3:
            print(f"{n:2d} | {nstr(fabs(x)**(mpf(1)/n),8)} | {nstr(fabs(H[n])**(mpf(1)/n),8)} | "
                  f"{nstr(En*(2*rho)**n,6)} | {ok} | {nstr(Nn,8)} | {nstr(law1,8)} | {nstr(res,6)}")
        assert ok
    msq = sum(Nn ** 2 for (n, Nn, _) in rows if n >= 20) / len([1 for (n, _, _) in rows if n >= 20])
    msq_law = sum(l ** 2 for (n, _, l) in rows if n >= 20) / len([1 for (n, _, _) in rows if n >= 20])
    pred = 2 * rho / pi if A == B else rho / pi
    print(f"mean N_n^2 (20<=n<=40) = {nstr(msq,6)}; mean law1^2 = {nstr(msq_law,6)}; Lemma 6 limit = {nstr(pred,6)}")

    # Lemma 1: t-form vs Laplace form at several eta
    for eta in [mpf("0.3"), mpf("-0.3"), mpf("0.7")]:
        tform = quad(lambda t: sqrt(A * A * t ** (2 * eta) + B * B * (1 - t) ** (2 * eta)), [0, mpf(1)/2, 1]) / r
        K = lambda u: (1 + exp(-u)) ** (-2 - eta)
        uform = quad(lambda u: exp(-u) * Hfun(A, B, eta * u) * K(u), [0, 5, 20, inf])
        print(f"Lemma 1, eta={nstr(eta,3)}: t-form {nstr(tform,20)}  Laplace {nstr(uform,20)}")
