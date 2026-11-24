# Referee check of thm:P2 (exact length at every angle). Own code, no project imports.
# L_theta(n) = int_0^1 |u + n x^{n-1} v| dx, u.v = c = cos(theta).
from mpmath import mp, mpf, quad, sqrt, cos, sin, pi, hyp2f1, log, exp, inf, nstr, fabs
import sys

mp.dps = 60

def legendre_list(M, c):
    P = [mpf(1), c]
    for k in range(1, M + 1):
        P.append(((2 * k + 1) * c * P[k] - k * P[k - 1]) / (k + 1))
    return P

def Lam_eps(eps, th):
    z = sin(th / 2) ** 2
    return 2 * z * hyp2f1(-eps, 1 + eps, 2, z)

def formula(n, th):
    n = mpf(n); c = cos(th); eps = 1 / (n - 1)
    head = 1 + c + pi * eps / sin(pi * eps) * Lam_eps(eps, th) * n ** (-eps)
    # number of terms so that n^-m < 10^-(dps+10)
    M = int((mp.dps + 15) * log(10) / log(n)) + 5
    P = legendre_list(M + 2, c)
    s = mpf(0)
    for m in range(1, M + 1):
        Lm = (P[m - 1] - P[m + 1]) / (2 * m + 1)
        s += (-1) ** m * eps / (m - eps) * Lm * n ** (-m)
    return head + s

def direct(n, th):
    # tau = -log x, dx = e^{-tau} dtau, slope s = n e^{-(n-1) tau}
    n = mpf(n); c = cos(th)
    def f(t):
        s = n * exp(-(n - 1) * t)
        return sqrt(1 + 2 * c * s + s * s) * exp(-t)
    ts = log(n) / (n - 1)
    w = 1 / (n - 1)
    pts = sorted(set([mpf(0)] + [ts + k * w for k in range(-40, 41) if ts + k * w > 0] + [ts + 40 * w + 2 ** j for j in range(0, 8)]))
    return quad(f, pts) + quad(f, [pts[-1], inf])

def mellin_lhs(eps, th):
    c = cos(th)
    f = lambda s: (sqrt(1 + 2 * c * s + s * s) - s - c) * s ** (eps - 1)
    return quad(f, [0, mpf(1) / 4, 1, 4, 64, inf])

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s); sys.stdout.flush()

pr("== Mellin identity: int (r-s-c)s^{eps-1} = pi/sin(pi eps) Lambda_eps ==")
worst = 0
for th in [mpf('0.3'), pi / 3, pi / 2, 2 * pi / 3, mpf('2.9')]:
    for eps in [mpf('0.1'), mpf('0.37'), mpf('0.5'), mpf('0.83')]:
        a = mellin_lhs(eps, th); b = pi / sin(pi * eps) * Lam_eps(eps, th)
        d = fabs(a - b) / fabs(b); worst = max(worst, d)
        pr(" th=%s eps=%s rel diff %s" % (nstr(th, 6), nstr(eps, 3), nstr(d, 3)))
pr(" worst", nstr(worst, 3))

# Lambda_eps alternative form: integral of Legendre P_eps from cos th to 1
from mpmath import legenp
pr("== Lambda_eps = int_c^1 P_eps(t)dt ==")
for th in [pi / 3, mpf('2.9')]:
    for eps in [mpf('0.37'), mpf('1.6')]:
        a = quad(lambda t: legenp(eps, 0, t, type=2), [cos(th), 1]); b = Lam_eps(eps, th)
        pr(" th=%s eps=%s diff %s" % (nstr(th, 5), nstr(eps, 3), nstr(fabs(a - b), 3)))

pr("== exact length: direct quadrature vs formula ==")
worst = 0
for th in [mpf('0.3'), pi / 3, pi / 2, 2 * pi / 3, 3 * pi / 4, mpf('2.9'), mpf('3.1')]:
    for n in ['1.3', '1.7', '2.5', '3', '7.5', '100', '10000', '1e8']:
        a = direct(n, th); b = formula(n, th)
        d = fabs(a - b); worst = max(worst, d)
        pr(" th=%s n=%s L=%s |diff|=%s" % (nstr(th, 5), n, nstr(a, 20), nstr(d, 3)))
pr(" worst |diff|", nstr(worst, 3))

pr("== removable singularities n=2 (eps=1), n=1.5 (eps=2) ==")
mp.dps = 90
for th in [pi / 3, mpf('2.9')]:
    for n0 in [2, mpf(3) / 2]:
        a = direct(n0, th)
        b1 = formula(mpf(n0) + mpf('1e-30'), th); b2 = formula(mpf(n0) - mpf('1e-30'), th)
        pr(" th=%s n=%s direct-formula(n0+1e-30)=%s , (n0-1e-30)=%s" % (nstr(th, 5), nstr(n0, 3), nstr(a - b1, 3), nstr(a - b2, 3)))
mp.dps = 60

pr("== theta=pi/2: check against the Ahmed form quoted in notes (not Ahmed's text) ==")
from mpmath import gamma, binomial
for n in [3, 10, 100]:
    n = mpf(n); eps = 1 / (n - 1)
    R = gamma(1 + eps / 2) * gamma((1 - eps) / 2) / (sqrt(pi) * (1 + eps))
    s = sum(binomial(mpf(1) / 2, k) / (n ** (2 * k - 1) * ((2 * k - 1) * (n - 1) - 1)) for k in range(1, 200))
    pr(" n=%s formula - [1+R n^-eps - sum] = %s" % (nstr(n, 4), nstr(formula(n, pi / 2) - (1 + R * n ** (-eps) - s), 3)))

pr("== Lambda_1 = sin^2/2 ; |Lambda_m| <= 1-c spot check ==")
for th in [pi / 3, mpf('2.9')]:
    c = cos(th); P = legendre_list(60, c)
    L1 = (P[0] - P[2]) / 3
    mx = max(fabs((P[m - 1] - P[m + 1]) / (2 * m + 1)) for m in range(1, 59))
    pr(" th=%s Lambda1-sin^2/2=%s max|Lambda_m|=%s 1-c=%s" % (nstr(th, 4), nstr(L1 - sin(th) ** 2 / 2, 3), nstr(mx, 6), nstr(1 - c, 6)))

open(__file__.replace('.py', '.out'), 'w').write("\n".join(out) + "\n")
