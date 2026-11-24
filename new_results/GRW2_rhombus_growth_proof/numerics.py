"""GRW2 numerics: Taylor coefficients I_n of I(eta) = P/(4r(1+eta)) for rhombi |x/A|^p+|y/B|^p=1.

I_n = sum_{j=0}^n H_{n-j} M_{n,j}           (Lemma 2 of proof.md)
H_m  = [z^m] (sqrt(A^2+B^2 e^{-2z}) + sqrt(B^2+A^2 e^{-2z}))/r      (power-series sqrt, high precision)
M_{n,0} = n! eta_D(n) exactly;  M_{n,j} (j>=1) = int_0^oo e^{-u} u^{n-j} k_j(u) du by tanh-sinh quadrature,
k_j(u) = (1+e^{-u})^{-2} (-log(1+e^{-u}))^j / j!.

Checks printed:
 (1) |I_n/n!|^{1/n} vs 1/rho;
 (2) I_n/n! against the one- and two-term laws of Theorem B;
 (3) E_n = I_n/n! - H_n against the bound (2 sigma)^{-n} (Lemma 4);
 (4) INDEPENDENT check of the coefficients: optimally truncated sum_n I_n eta^n vs I(eta) from
     direct quadrature of the perimeter arc-length integral, at eta = +-0.04.
Usage: python numerics.py NMAX
"""
import sys
import mpmath as mp

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp.mp.dps = 45
CASES = [(1, 1), (1, 2), (1, 5), (2, 3)]


def sqrt_series(a, N):
    s = [mp.sqrt(a[0])] + [mp.mpf(0)] * N
    for m in range(1, N + 1):
        acc = a[m] - mp.fsum(s[k] * s[m - k] for k in range(1, m))
        s[m] = acc / (2 * s[0])
    return s


def H_coeffs(A, B, N):
    A, B = mp.mpf(A), mp.mpf(B)
    r = mp.sqrt(A**2 + B**2)
    e = [(-2) ** m / mp.factorial(m) for m in range(N + 1)]  # e^{-2z}
    a1 = [B**2 * c for c in e]; a1[0] += A**2
    a2 = [A**2 * c for c in e]; a2[0] += B**2
    s1, s2 = sqrt_series(a1, N), sqrt_series(a2, N)
    return [(x + y) / r for x, y in zip(s1, s2)]


def M(n, j):
    if j == 0:
        if n == 0:
            return mp.mpf(1) / 2
        if n == 1:
            return mp.log(2)
        return mp.factorial(n) * mp.altzeta(n)
    k = n - j
    f = lambda u: mp.exp(-u) * u**k * (-mp.log1p(mp.exp(-u)))**j / (1 + mp.exp(-u))**2
    peak = mp.mpf(k) / (j + 1)
    pts = sorted(set([mp.mpf(0), peak / 2, peak, 2 * peak + 1, 4 * peak + 8, 8 * peak + 40, mp.inf]))
    return mp.quad(f, pts) / mp.factorial(j)


def I_direct(A, B, eta):
    """I(eta) = P/(4r(1+eta)) from the arc-length integral, no series used."""
    A, B = mp.mpf(A), mp.mpf(B)
    r = mp.sqrt(A**2 + B**2)
    f = lambda t: mp.sqrt(A**2 * t**(2 * eta) + B**2 * (1 - t)**(2 * eta))
    return mp.quad(f, [0, mp.mpf(1) / 2, 1]) / r


def laws(A, B, n):
    A, B = mp.mpf(A), mp.mpf(B)
    r = mp.sqrt(A**2 + B**2)
    L = mp.log(B / A)
    z0 = mp.mpc(L, mp.pi / 2)
    z1 = mp.mpc(-L, mp.pi / 2)
    pts = [(z0, A), (mp.conj(z0), A), (z1, B), (mp.conj(z1), B)]
    one = two = mp.mpc(0)
    for zj, a in pts:
        c = a / r * mp.sqrt(-2 * zj)  # principal branch, Lemma 3
        base = -c * zj**(-n) * mp.mpf(n)**(-1.5) / (2 * mp.sqrt(mp.pi))
        one += base
        two += base * (1 + mp.mpf(3) / (8 * n) - 3 * zj / (4 * n))
    return mp.re(one), mp.re(two)


def bracket_form(A, B, n):
    A, B = mp.mpf(A), mp.mpf(B)
    r = mp.sqrt(A**2 + B**2)
    L = mp.log(B / A)
    rho = mp.sqrt(L**2 + mp.pi**2 / 4)
    th = mp.atan2(mp.pi / 2, L)
    br = B * (-1)**n * mp.cos((n - mp.mpf(1) / 2) * th) - A * mp.sin((n - mp.mpf(1) / 2) * th)
    return -mp.sqrt(2 * rho) / (r * mp.sqrt(mp.pi)) * rho**(-n) * mp.mpf(n)**(-1.5) * br, br


def run(A, B):
    L = mp.log(mp.mpf(B) / A)
    rho = mp.sqrt(L**2 + mp.pi**2 / 4)
    r = mp.sqrt(mp.mpf(A)**2 + B**2)
    print(f"\n=== (A,B)=({A},{B})  rho={mp.nstr(rho,15)}  1/rho={mp.nstr(1/rho,12)}  theta0/pi={mp.nstr(mp.atan2(mp.pi/2,L)/mp.pi,10)}")
    H = H_coeffs(A, B, NMAX)
    I = []
    for n in range(NMAX + 1):
        I.append(mp.fsum(H[n - j] * M(n, j) for j in range(n + 1)))
    # (1)-(3)
    print(" n | |I_n/n!|^(1/n) | I_n/n! * rho^n n^1.5 | law1 (same norm) | (I/n!-law2)/env*n^2 | E_n*(1.5rho)^n | bracket")
    worst2 = mp.mpf(0)
    for n in range(1, NMAX + 1):
        b = I[n] / mp.factorial(n)
        l1, l2 = laws(A, B, n)
        lb, br = bracket_form(A, B, n)
        assert abs(lb - l1) < mp.mpf(10)**(-35) * rho**(-n)  # bracket form == sum over 4 points
        env = rho**(-n) * mp.mpf(n)**(-1.5)
        E = b - H[n]
        d2 = (b - l2) / env * n**2
        if n >= 10:
            worst2 = max(worst2, abs(d2))
        if n % 5 == 0 or n <= 4:
            print(f"{n:3d} | {mp.nstr(abs(b)**(mp.mpf(1)/n),8):>10} | {mp.nstr(b/env,10):>14} | {mp.nstr(l1/env,10):>14} | {mp.nstr(d2,6):>10} | {mp.nstr(E*(mp.mpf(1.5)*rho)**n,5):>10} | {mp.nstr(br,6)}")
    print(f" max_(10<=n<={NMAX}) |I_n/n! - law2| rho^n n^(3.5) = {mp.nstr(worst2,6)}")
    # mean square of normalized coefficients vs r^2/2 * (2 rho/pi)/r^2 = rho/pi  (Lemma 6 prediction for L!=0)
    ms = mp.fsum((I[n] / mp.factorial(n) * rho**n * mp.mpf(n)**1.5)**2 for n in range(NMAX // 2, NMAX + 1)) / (NMAX - NMAX // 2 + 1)
    pred = rho / mp.pi if A != B else 2 * rho / mp.pi
    print(f" mean of (I_n/n! rho^n n^1.5)^2 over n in [{NMAX//2},{NMAX}] = {mp.nstr(ms,6)}  (limit predicted {mp.nstr(pred,6)})")
    # (4) independent check
    for eta in (mp.mpf('0.04'), mp.mpf('-0.04')):
        ref = I_direct(A, B, eta)
        best = None
        partial = mp.mpf(0)
        for n in range(NMAX + 1):
            partial += I[n] * eta**n
            err = abs(partial - ref)
            if best is None or err < best[1]:
                best = (n, err)
        smallest = min(abs(I[n] * eta**n) for n in range(NMAX + 1))
        print(f" eta={mp.nstr(eta,3)}: I(eta) direct = {mp.nstr(ref,25)}; best truncation N={best[0]}, |error|={mp.nstr(best[1],3)}; smallest term {mp.nstr(smallest,3)}; e^(-rho/|eta|)={mp.nstr(mp.exp(-rho/abs(eta)),3)}")
    return I


if __name__ == "__main__":
    for A, B in CASES:
        run(A, B)
