# Referee check of thm:fn (length of the graph of f^n). Own code.
from mpmath import mp, mpf, quad, sqrt, log, exp, sin, cos, pi, e, zeta, inf, nstr, fabs, expm1
import sys
mp.dps = 60
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s); sys.stdout.flush()

# ---- random slope cumulants ----
gp = lambda s: 1 - s / sqrt(1 + s * s)
mom = [quad(lambda s: gp(s) * (-log(s)) ** k, [0, mpf(1) / 8, 1, 8, inf]) for k in range(5)]
k1 = mom[1]; m2 = mom[2] - k1 ** 2
c3 = mom[3] - 3 * mom[2] * k1 + 2 * k1 ** 3
mu = lambda k: quad(lambda s: gp(s) * (-log(s) - k1) ** k, [0, mpf(1) / 8, 1, 8, inf])
cm2, cm3, cm4 = mu(2), mu(3), mu(4)
kk4 = cm4 - 3 * cm2 ** 2
pr("mass of g'", nstr(mom[0] - 1, 3))
pr("k1 - (1-log2) =", nstr(k1 - (1 - log(2)), 3))
pr("k2 - (1+pi^2/6) =", nstr(cm2 - (1 + pi ** 2 / 6), 3))
pr("k3 - (2-1.5 zeta3) =", nstr(cm3 - (2 - mpf(3) / 2 * zeta(3)), 3))
pr("k4 - (6+pi^4/15) =", nstr(kk4 - (6 + pi ** 4 / 15), 3))
pr("geometric mean of S - 2/e =", nstr(exp(-k1) - 2 / e, 3))
K1, K2, K3 = 1 - log(2), 1 + pi ** 2 / 6, 2 - mpf(3) / 2 * zeta(3)

# ---- test functions: (f, f', a, A=f''(1), B=f'''(1)) ----
E1 = e - 1
funcs = {
 'x': (lambda x: x, lambda x: mpf(1), mpf(1), mpf(0), mpf(0)),
 '(e^x-1)/(e-1)': (lambda x: expm1(x) / E1, lambda x: exp(x) / E1, e / E1, e / E1, e / E1),
 '2x/(1+x)': (lambda x: 2 * x / (1 + x), lambda x: 2 / (1 + x) ** 2, mpf(1) / 2, -mpf(1) / 2, mpf(3) / 4),
 'sin(pi x/4)/sin(pi/4)': (lambda x: sin(pi * x / 4) / sin(pi / 4), lambda x: pi / 4 * cos(pi * x / 4) / sin(pi / 4), pi / 4, -(pi / 4) ** 2, -(pi / 4) ** 3),
 '(x+x^3)/2': (lambda x: (x + x ** 3) / 2, lambda x: (1 + 3 * x ** 2) / 2, mpf(2), mpf(3), mpf(3)),
}

def defect(F, Fp, a, n, scale=None):
    # D = int_0^1 g(y') dx, g(s) = 1 - 1/(s + sqrt(1+s^2)); integrate in tau = 1-x
    n = mpf(n)
    def h(t):
        x = 1 - t
        fx = F(x)
        if fx <= 0: return mpf(0)
        s = n * exp((n - 1) * log(fx)) * Fp(x)
        return 1 - 1 / (s + sqrt(1 + s * s))
    w = scale if scale is not None else 1 / (n * a)
    pts = [mpf(0)]
    k = mpf(1) / 4
    while k * w < 1:
        pts.append(k * w); k *= mpf(3) / 2
    pts.append(mpf(1))
    return quad(h, pts)

def J(a, A, B, n):
    n = mpf(n); eps = 1 / (n - 1); z = mpf(1) / 2; c = 0
    Xi = log(n * a) + K1
    m1 = Xi; m2 = Xi ** 2 + K2; m3 = Xi ** 3 + 3 * K2 * Xi + K3
    J1 = 2 * z / a * m1
    J2 = z / a ** 3 * (A * (m2 - 2 * m1) - a ** 2 * m2 + 2 * a * (1 - z))
    J3 = z / (3 * a ** 5) * (a ** 4 * m3 - (3 * A * a ** 2 + B * a) * (m3 - 3 * m2) + 3 * A ** 2 * (m3 - 4 * m2 + 2 * m1)) - 4 * A * z * (1 - z) / a ** 4 - c * z * (1 - z) / a ** 3
    return eps, J1, J2, J3

pr("== Thm fn: D_n=2-L_n vs (log n+log a+1-log2)/(a(n-1)) ; and eps J1+eps^2J2+eps^3J3 (theta=pi/2) ==")
for name, (F, Fp, a, A, B) in funcs.items():
    pr("f =", name, " a=", nstr(a, 8), " A=", nstr(A, 8), " B=", nstr(B, 8))
    for n in ['1e3', '1e4', '1e5', '1e6', '1e8', '1e10']:
        D = defect(F, Fp, a, n)
        nn = mpf(n)
        main = (log(nn) + log(a) + 1 - log(2)) / (a * (nn - 1))
        eps, J1, J2, J3 = J(a, A, B, n)
        r1 = (D - main) / (log(nn) ** 2 / nn ** 2)
        r3 = (D - eps * J1 - eps ** 2 * J2 - eps ** 3 * J3) / (eps ** 4 * log(nn) ** 4)
        pr("  n=%s D=%s  (D-main)/(log^2n/n^2)=%s  (D-J1..J3)/(eps^4 log^4 n)=%s" % (n, nstr(D, 25), nstr(r1, 8), nstr(r3, 6)))
    pr("  predicted limit of first ratio: -(a^2-A)/(2a^3) =", nstr(-(a ** 2 - A) / (2 * a ** 3), 8))

# accuracy check of quadrature: two precisions
mp.dps = 80
F, Fp, a, A, B = funcs['(e^x-1)/(e-1)']
D80 = defect(F, Fp, a, '1e6')
mp.dps = 60
D60 = defect(F, Fp, a, '1e6')
pr("quadrature self-check (e^x-1)/(e-1), n=1e6: |D60-D80| =", nstr(fabs(D60 - D80), 3))

# small n: is the O(log^2 n/n^2) remainder already tame at moderate n?  (just informative)
pr("== hypothesis probe: f in C^1 only away from 1 (|x-1/2|^{3/2}-type kink), C^2 near 1 ==")
# f(x) = (x + 0.2*|x-1/2|^{1.5}*sign)/norm ; nondecreasing, C^1, f(1)=1
def Fk(x):
    t = x - mpf(1) / 2
    return (x + mpf('0.2') * (fabs(t) ** mpf('1.5')) * (1 if t >= 0 else -1) + mpf('0.2') * mpf('0.5') ** mpf('1.5')) / (1 + mpf('0.4') * mpf('0.5') ** mpf('1.5'))
def Fkp(x):
    t = x - mpf(1) / 2
    return (1 + mpf('0.3') * fabs(t) ** mpf('0.5')) / (1 + mpf('0.4') * mpf('0.5') ** mpf('1.5'))
nrm = 1 + mpf('0.4') * mpf('0.5') ** mpf('1.5')
ak = (1 + mpf('0.3') * mpf('0.5') ** mpf('0.5')) / nrm
pr(" f(0)=", nstr(Fk(mpf(0)), 5), " f(1)=", nstr(Fk(mpf(1)), 5))
for n in ['1e4', '1e6', '1e8']:
    nn = mpf(n)
    D = defect(Fk, Fkp, ak, n)
    main = (log(nn) + log(ak) + 1 - log(2)) / (ak * (nn - 1))
    pr("  n=%s (D-main)/(log^2n/n^2)=%s" % (n, nstr((D - main) / (log(nn) ** 2 / nn ** 2), 8)))

# ---- contact law at theta=pi/2: f = sin(pi x/2), m=2, beta=pi^2/8 ----
pr("== contact law f=sin(pi x/2): D vs (log n/(m beta n))^(1/m) and the 1/Lambda-corrected form ==")
F = lambda x: sin(pi * x / 2); Fp = lambda x: pi / 2 * cos(pi * x / 2)
m = 2; beta = pi ** 2 / 8
for n in ['1e4', '1e6', '1e8', '1e12', '1e16']:
    nn = mpf(n); N = (beta * nn) ** (mpf(1) / m); Lam = log(m * N)
    D = defect(F, Fp, mpf(1), n, scale=1 / N)
    lead = (log(nn) / (m * beta * nn)) ** (mpf(1) / m)
    corr = Lam ** (mpf(1) / m) / N * (1 + ((m - 1) * log(Lam) + m * K1) / (m ** 2 * Lam))
    pr("  n=%s D/lead=%s  D/corrected=%s" % (n, nstr(D / lead, 8), nstr(D / corr, 8)))

open(__file__.replace('.py', '.out'), 'w').write("\n".join(out) + "\n")
