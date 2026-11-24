"""pade_borel_check.py -- the divergent diamond-end series is Borel summable (numerical check).

Paper: "Investigating zeta(3)", Section XI-A ("The diamond end: radius zero"), the sentence
"a Pade-Borel resummation built from f_0, ..., f_180 reproduces direct quadrature of I(eta)
at eta = 0.1, 0.2, 0.3 and -0.2 to 44, 31, 25 and 30 digits".

For the square, I(eta) = P(p) / (4 sqrt2 (1+eta)) = int_0^1 (t(1-t))^(eta/2) sqrt(cosh(eta l)) dt,
l = log(t/(1-t)), eta = 1/p - 1, and f_n = [eta^n] I.  The f_n grow like n! (2/pi)^n, so the series
diverges (Theorem "radius 0").  This script

  1. computes f_0 .. f_180 at 450 digits from the Gamma-function generating function
     (log Gamma(1+x) = -gamma x + sum_n (-1)^n zeta(n) x^n / n; Euler's gamma cancels),
  2. builds the [N/N] Pade approximant (N = 80, 85, 90) of the Borel transform of the subtracted
     coefficients f_n - sqrt2 (-1)^n, so that I(eta) = sqrt2/(1+eta) + int_0^inf e^-s PB(eta s) ds,
  3. compares the Laplace integral with a fresh direct quadrature of I(eta) at 110 digits.

Origin: research job D11b (Kimi K3), steps 29-30; repackaged here with no file I/O and a fresh
direct quadrature.  Run:  python pade_borel_check.py      (a few minutes; prints digits per eta)
"""
import time
from math import factorial
from fractions import Fraction
from mpmath import mp, mpf, zeta, binomial, sqrt, quad, exp, log, cosh, nstr

T0 = time.time()
NMAX = 180


def coefficients(nmax=NMAX, dps=450):
    """f_0..f_nmax of I(eta) for the square, from the Gamma generating function (exact recursion)."""
    mp.dps = dps
    Zv = [mpf(0)] * 2 + [zeta(n) for n in range(2, nmax + 3)]
    A = {}
    for j in range(0, nmax + 1, 2):
        arr = [mpf(0)] * (nmax + 1)
        for m in range(0, nmax + 1 - j):
            n = j + m
            if n >= 2:
                arr[m] = 2 * ((-1) ** n) * Zv[n] / n * binomial(n, j)
        A[j] = arr

    def conv(a, b, N):
        c = [mpf(0)] * (N + 1)
        for i in range(N + 1):
            if a[i] == 0:
                continue
            for j in range(N + 1 - i):
                if b[j] != 0:
                    c[i + j] += a[i] * b[j]
        return c

    H0 = [mpf(1)] + [mpf(0)] * nmax
    for m in range(1, nmax + 1):
        H0[m] = sum(i * A[0][i] * H0[m - i] for i in range(2, m + 1)) / m
    Hs = {0: H0}
    for j in range(2, nmax + 1, 2):
        acc = [mpf(0)] * (nmax + 1)
        for i in range(2, j + 1, 2):
            sc = conv(A[i], Hs[j - i], nmax)
            for m in range(nmax + 1):
                acc[m] += i * sc[m]
        Hs[j] = [x / j for x in acc]
    # Taylor coefficients of sqrt(cosh z) = sum cc_k z^(2k)
    ccf = [Fraction(1)]
    for k in range(1, nmax // 2 + 1):
        ccf.append((Fraction(1, factorial(2 * k)) - sum(ccf[i] * ccf[k - i] for i in range(1, k))) / 2)
    cc = [mpf(x.numerator) / mpf(x.denominator) for x in ccf]
    main = [mpf(0)] * (nmax + 1)
    for k in range(0, nmax // 2 + 1):
        for m in range(0, nmax + 1 - 2 * k):
            main[2 * k + m] += cc[k] * Hs[2 * k][m] * factorial(2 * k) / mpf(2) ** m
    logK = [mpf(0)] * (nmax + 1)
    for n in range(2, nmax + 1):
        logK[n] = -((-1) ** n) * Zv[n] / n
    Ks = [mpf(1)] + [mpf(0)] * nmax
    for m in range(1, nmax + 1):
        Ks[m] = sum(j * logK[j] * Ks[m - j] for j in range(2, m + 1)) / m
    inv1p = [mpf((-1) ** m) for m in range(nmax + 1)]
    return conv(conv(Ks, inv1p, nmax), main, nmax)


def pade(a, L, M):
    """[L/M] Pade approximant of sum a_n t^n (Hankel system)."""
    A = mp.matrix(M, M)
    rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            k = L + i - j
            A[i, j] = a[k] if k >= 0 else mpf(0)
        rhs[i] = -a[L + 1 + i]
    q = mp.lu_solve(A, rhs)
    den = [mpf(1)] + [q[j] for j in range(M)]
    num = [sum(den[j] * a[n - j] for j in range(0, min(M, n) + 1)) for n in range(L + 1)]
    return num, den


def peval(num, den, t):
    nv = mp.polyval(num[::-1], t)
    dv = mp.polyval(den[::-1], t)
    return nv / dv


def laplace(num, den, eta, dps=100):
    old = mp.dps
    mp.dps = dps
    try:
        pts = [0, 1, 4, 16, 40, mp.inf] if eta > 0 else list(range(0, 84, 4)) + [84]
        return quad(lambda s: exp(-s) * peval(num, den, eta * s), pts)
    finally:
        mp.dps = old


def direct(eta, dps=110):
    old = mp.dps
    mp.dps = dps
    try:
        f = lambda t: (t * (1 - t)) ** (eta / 2) * sqrt(cosh(eta * log(t / (1 - t))))
        return quad(f, [0, mpf(1) / 4, mpf(1) / 2, mpf(3) / 4, 1])
    finally:
        mp.dps = old


if __name__ == "__main__":
    f = coefficients()
    print(f"f_0..f_{NMAX} at 450 digits: {time.time() - T0:.0f}s;  f_1 = {nstr(f[1], 5)}, "
          f"f_2 = {nstr(f[2], 20)} (pi^2/24 = {nstr(mp.pi**2 / 24, 20)})", flush=True)
    mp.dps = 400
    s2 = sqrt(2)
    borel = [(f[n] - s2 * (-1) ** n) / mpf(factorial(n)) for n in range(len(f))]
    P = {N: pade(borel[:2 * N + 1], N, N) for N in (80, 85, 90)}
    print(f"Pade approximants built: {time.time() - T0:.0f}s", flush=True)
    for es in ("0.1", "0.2", "0.3", "-0.2"):
        eta = mpf(es)
        ref = direct(eta)
        vals = {N: s2 / (1 + eta) + laplace(*P[N], eta) for N in P}
        dig = -mp.log10(abs(vals[90] - ref) / abs(ref))
        spread = max(abs(vals[90] - vals[85]), abs(vals[85] - vals[80]))
        print(f"eta = {es:>5}:  Pade-Borel [90/90] = {nstr(vals[90], 50)}\n"
              f"               direct quadrature = {nstr(ref, 50)}\n"
              f"               agreement {nstr(dig, 4)} digits   ([80/80]..[90/90] spread {nstr(spread, 3)})",
              flush=True)
    print(f"total {time.time() - T0:.0f}s")
