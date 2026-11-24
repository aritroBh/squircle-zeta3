"""High-precision numerics for M_k (k = 0..3), three independent methods.

 Method A  (Euler variable v):  M_k = int_b^a sum_{ij} R_ij(v) log^i t log^j (1-t) dv  with the exact rational
           R_ij from mk_integrands.py; tanh-sinh quadrature.  The R_ij vanish where the logs blow up, so
           the integrand is bounded and smooth except for log-power factors at the endpoints.
 Method B  (series under the integral, independent parametrisation, as in research/05_final_verify):
           x = s^u, y = (1-s)^u; P/8 = 2^{-u} + u J(u), J(u) = int_0^{1/2} s^{u-1}(sqrt(1+r^{2-2u}) - 1) ds,
           r = s/(1-s);  M_n = J_n + (-log 2)^{n+1}/(n+1)!, with J_n from the Taylor coefficients in u
           (computed symbolically here).
 Method C  (black box): P(p) by direct quadrature at large p, subtract 8(1+M0/p+M1/p^2+M2/p^3) using the
           proved closed forms M0, M1 and the method-A value of M2, Neville-extrapolate p^4(...) -> 8 M3.

Usage:  python mk_numerics.py A <dps> | B <dps> | C <dps>
Writes ../outputs/M{k}_method{X}_{dps}.txt
"""
import os, sys, time
import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "outputs")
method = sys.argv[1]
DPS = int(sys.argv[2])
KMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 3
mp.mp.dps = DPS + 30
T0 = time.time()

r2 = mp.sqrt(2); aN = 1 + r2; bN = r2 - 1; lam = mp.log(2); LN = mp.log(aN)
M0_cf = lam - r2 * LN
M1_cf = (-2 * mp.pi ** 2 + 19 * r2 * mp.pi ** 2 + 24 * lam ** 2 - 48 * LN ** 2 - 120 * r2 * LN ** 2
         + 96 * r2 * lam * LN - 240 * r2 * mp.polylog(2, bN)) / 48


def save(name, val):
    with open(os.path.join(OUT, f"{name}_method{method}_{DPS}.txt"), "w") as fh:
        fh.write(mp.nstr(val, DPS) + "\n")


def report(k, val, err=None):
    with mp.workdps(DPS):
        print(f"M{k} = {mp.nstr(val, 60)}...   (quad err est {mp.nstr(err, 3) if err is not None else '-'};"
              f" {time.time() - T0:.0f}s)", flush=True)


if method == "A":
    from mk_integrands import v_integrand, v as vs
    for k in range(0, KMAX + 1):
        Rk = v_integrand(k)
        fs = [(i, j, sp.lambdify(vs, R, modules="mpmath")) for (i, j), R in Rk.items()]

        def f(x, fs=fs):
            l1 = mp.log((x - bN) * (x + aN) / (4 * x))
            l2 = mp.log((aN - x) * (x + bN) / (4 * x))
            return mp.fsum(R(x) * l1 ** i * l2 ** j for i, j, R in fs)

        val, err = mp.quad(f, [bN, mp.mpf(1) / 2, 1, 2, aN], error=True, maxdegree=14)
        report(k, val, err)
        save(f"M{k}", val)
        if k == 0:
            print("   M0 - (log2 - sqrt2 L) =", mp.nstr(val - M0_cf, 3))
        if k == 1:
            print("   M1 - closed form      =", mp.nstr(val - M1_cf, 3))

elif method == "B":
    u, w, l, lm = sp.symbols('u w l lm')
    # sqrt(1 + w e^{-2 u l}) - 1 and s^{u-1} = s^{-1} e^{u lm};  Taylor in u
    expr = (sp.sqrt(1 + w * sp.exp(-2 * u * l)) - 1) * sp.exp(u * lm)
    ser = sp.series(expr, u, 0, KMAX + 1).removeO()
    coeffs = [sp.simplify(ser.coeff(u, n)) for n in range(KMAX + 1)]
    fns = [sp.lambdify((w, l, lm), c, modules="mpmath") for c in coeffs]
    for n in range(KMAX + 1):
        def f(s, F=fns[n]):
            r = s / (1 - s)
            return F(r * r, mp.log(r), mp.log(s)) / s
        Jn, err = mp.quad(f, [0, mp.mpf(1) / 256, mp.mpf(1) / 32, mp.mpf(1) / 8, mp.mpf(1) / 4, mp.mpf(1) / 2],
                          error=True, maxdegree=14)
        Mn = Jn + (-lam) ** (n + 1) / mp.factorial(n + 1)
        report(n, Mn, err)
        save(f"M{n}", Mn)

elif method == "C":
    # need M2 from method A at higher precision
    M2 = mp.mpf(open(os.path.join(OUT, "M2_methodA_300.txt")).read().strip())

    def F(eps):
        def phi(t):
            a = t ** (eps - 1); b = (1 - t) ** (eps - 1)
            z = (b / a) ** 2
            return b - b * b / (a * (mp.sqrt(1 + z) + 1))
        return 2 * mp.quad(phi, [0, mp.mpf(1) / 64, mp.mpf(1) / 8, mp.mpf(1) / 2])

    ps = [mp.mpf(400) * 2 ** j for j in range(10)]
    ys = []
    for p in ps:
        uu = 1 / p
        Mu = -F(uu) / 2                     # M(u) = -F(u)/2  (P/8 = 1 + u M(u))
        ys.append((Mu - M0_cf - M1_cf * uu - M2 * uu ** 2) / uu ** 3)    # -> M3 + M4 u + ...
        print(f"  p = {int(p)}:  (M(u)-M0-M1u-M2u^2)/u^3 = {mp.nstr(ys[-1], 30)}  ({time.time()-T0:.0f}s)", flush=True)
    # Neville extrapolation to u = 0
    xs = [1 / p for p in ps]
    n = len(xs)
    P = list(ys)
    for m in range(1, n):
        for i in range(n - m):
            P[i] = ((0 - xs[i + m]) * P[i] + (xs[i] - 0) * P[i + 1]) / (xs[i] - xs[i + m])
    print("Neville extrapolation M3 =", mp.nstr(P[0], 40))
    save("M3", P[0])
print(f"total {time.time() - T0:.0f}s")
