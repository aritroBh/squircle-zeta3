# Referee check of thm:S (superellipsoid |x|^p+|y|^p+|z|^p=1). Own code, no project imports.
from mpmath import (mp, mpf, quad, sqrt, log, exp, pi, gamma, zeta, nstr, fabs, inf, expm1,
                    polylog, atan, euler, diff)
import sys, time
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s); sys.stdout.flush()
    open(__file__.replace('.py', '.out'), 'w').write("\n".join(out) + "\n")

T0 = time.time()
# ---------- direct surface area from the definition ----------
def S_direct(p):
    p = mpf(p)
    # part of the first octant where z is the largest coordinate and y <= x; y = s x.
    def outer(s):
        xm = (2 + s ** p) ** (-1 / p)
        def inner(x):
            zp = 1 - x ** p * (1 + s ** p)
            z = zp ** (1 / p)
            return x * sqrt(1 + (x / z) ** (2 * p - 2) + (s * x / z) ** (2 * p - 2))
        return quad(inner, [0, xm / 2, xm * mpf('0.9'), xm])
    return 48 * quad(outer, [0, mpf(1) / 2, 1])

# ---------- squircle perimeter, cancellation-free ----------
def phi(t, e):
    A = exp((e - 1) * log(t)); B = exp((e - 1) * log(1 - t))
    return -2 * A * B / (sqrt(A * A + B * B) + A + B)
def I(e):
    return 2 * quad(lambda t: phi(t, e), [0, mpf(1) / 4, mpf(1) / 2])
def P(p):
    e = 1 / mpf(p); return 8 + 4 * e * I(e)

# ---------- vertex function Q ----------
def h(e, y):
    a = 2 - 2 * e
    f = lambda x: -expm1(-x ** a / y ** 2) * x ** (e - 1) * exp(-x)
    x0 = y ** (1 / a) if a > 0 else mpf(1)
    pts = sorted(set([q for q in [mpf(0), x0 / 4, x0, 4 * x0, mpf(1), mpf(8), mpf(40)] if q <= 40]))
    return quad(f, pts) + quad(f, [mpf(40), inf])
def H3(e):
    g = lambda y: h(e, y) ** 3
    return quad(g, [0, mpf(1) / 100, mpf(1) / 10, 1, 10, 100]) + quad(g, [100, inf])
def Q(e):
    e = mpf(e); return -H3(e) / (sqrt(pi) * gamma(1 + 2 * e))
def G(e):
    return gamma(1 + e) ** 2 / gamma(1 + 2 * e)
def S_formula(p):
    e = 1 / mpf(p)
    return 6 * G(e) * (P(p) - 4) - 8 * e ** 2 * Q(e)

lam = log(2); L = log(1 + sqrt(2)); b = sqrt(2) - 1; ac = atan(1 / sqrt(2))

mp.dps = 30
pr("== sanity: P(2)-2pi, P(1)-4sqrt2 ==", nstr(P(2) - 2 * pi, 3), nstr(P(1) - 4 * sqrt(2), 3))
pr("== exact formula vs direct surface integral ==")
for p in ['1', '2', '2.5', '3', '4', '10']:
    t = time.time()
    a = S_formula(p)
    if p == '1': d = 4 * sqrt(3)
    elif p == '2': d = 4 * pi
    else: d = S_direct(p)
    pr("  p=%s S_formula=%s S_ref=%s rel diff=%s (%.0fs)" % (p, nstr(a, 28), nstr(d, 28), nstr(fabs(a - d) / d, 3), time.time() - t))

# ---------- Q0 ----------
mp.dps = 40
Q0c = (-pi ** 2 / 8 + sqrt(2) / 16 * pi ** 2 + 3 * sqrt(2) / 2 * pi * ac - 3 * lam ** 2 - 9 * sqrt(2) * L * log(3)
       + (6 + 6 * sqrt(2)) * L ** 2 - 6 * sqrt(2) * polylog(2, b) + 9 * sqrt(2) / 4 * polylog(2, -b ** 4))
t = time.time(); Q0n = Q(0)
pr("== Q0: integral %s closed form %s diff %s (%.0fs)" % (nstr(Q0n, 35), nstr(Q0c, 35), nstr(Q0n - Q0c, 3), time.time() - t))

# ---------- Q1 = Q'(0) = -H3'(0) + 2*gamma_E*Q0, H3'(0) = (3/sqrt pi) int h0^2 dh ----------
def dh(y):
    f = lambda x: (-expm1(-x ** 2 / y ** 2) - 2 * (x ** 2 / y ** 2) * exp(-x ** 2 / y ** 2)) * log(x) * exp(-x) / x
    pts = sorted(set([q for q in [mpf(0), y / 4, y, 4 * y, mpf(1), mpf(8), mpf(40)] if q <= 40]))
    return quad(f, pts) + quad(f, [mpf(40), inf])
t = time.time()
g1 = lambda y: h(0, y) ** 2 * dh(y)
dH3 = 3 / sqrt(pi) * (quad(g1, [0, mpf(1) / 100, mpf(1) / 10, 1, 10, 100]) + quad(g1, [100, inf]))
Q1n = -dH3 + 2 * euler * Q0n
pr("== Q1 (derivative under the integral) %s ; printed 6.5043255737... diff %s (%.0fs)" % (nstr(Q1n, 30), nstr(Q1n - mpf('6.504325573741009758766107790505238481960823372240351285303526704017478'), 3), time.time() - t))
mp.dps = 25
dd = mpf('1e-4')
Q1fd = (Q(dd) - Q(-dd)) / (2 * dd)
pr("   crude central difference (delta=1e-4): %s (expected error ~1e-8)" % nstr(Q1fd, 15))
mp.dps = 40

# ---------- M0, M1, M2 from the perimeter (independent of the paper's closed forms) ----------
def Ider(k):
    return 2 * quad(lambda t: diff(lambda e: phi(t, e), 0, k), [0, mpf(1) / 4, mpf(1) / 2])
M0n = Ider(0) / 2; M1n = Ider(1) / 2; M2n = Ider(2) / 4
M0c = lam - sqrt(2) * L
M1c = (19 * sqrt(2) / 48 - mpf(1) / 24) * pi ** 2 + lam ** 2 / 2 + 2 * sqrt(2) * lam * L - (1 + 5 * sqrt(2) / 2) * L ** 2 - 5 * sqrt(2) * polylog(2, b)
pr("== M0 %s diff closed %s" % (nstr(M0n, 30), nstr(M0n - M0c, 3)))
pr("== M1 %s diff closed %s" % (nstr(M1n, 30), nstr(M1n - M1c, 3)))
pr("== M2 %s diff printed -0.02216191108726939942: %s" % (nstr(M2n, 30), nstr(M2n - mpf('-0.02216191108726939942'), 3)))

s1 = 48 * M0n
s2 = 48 * M1n - 4 * pi ** 2 - 8 * Q0n
s3 = 48 * M2n - 8 * pi ** 2 * M0n + 48 * zeta(3) - 8 * Q1n
pr("== s1 %s printed -26.5585583865" % nstr(s1, 25))
pr("== s2 %s printed -4.6892155184" % nstr(s2, 25))
pr("== s3 %s printed 48.2874324860" % nstr(s3, 25))
s2closed = (-5 * pi ** 2 + 37 * sqrt(2) / 2 * pi ** 2 - 12 * sqrt(2) * pi * ac + 48 * lam ** 2 - 96 * L ** 2 - 168 * sqrt(2) * L ** 2
            + 96 * sqrt(2) * lam * L + 72 * sqrt(2) * log(3) * L - 192 * sqrt(2) * polylog(2, b) - 18 * sqrt(2) * polylog(2, -b ** 4))
pr("   notes' boxed s2 closed form - s2: %s" % nstr(s2closed - s2, 3))

# ---------- end-to-end: S(p) at large p minus the three-term expansion ----------
mp.dps = 35
pr("== R4(p) = p^4 (S(p) - 24 - s1/p - s2/p^2 - s3/p^3)  (should approach s4 ~ -78.07, bounded) ==")
for p in [100, 300, 1000, 3000]:
    t = time.time()
    Sp = S_formula(p); pp = mpf(p)
    R4 = pp ** 4 * (Sp - 24 - s1 / pp - s2 / pp ** 2 - s3 / pp ** 3)
    pr("  p=%s S=%s R4=%s (%.0fs)" % (p, nstr(Sp, 30), nstr(R4, 10), time.time() - t))
pr("total time %.0fs" % (time.time() - T0))
