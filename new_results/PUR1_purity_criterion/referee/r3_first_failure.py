# Referee check of (e): first-failure formula via exact Taylor coefficients of F (Cauchy contour on the
# direct integral), the Euler-sum formula for gamma_{a,2}, and the classical weight-4 evaluation.
import mpmath as mp

out = []
P_ = lambda *a: (out.append(' '.join(str(t) for t in a)), print(*a))

# ---------- Euler sums by Cohen-Villegas-Zagier acceleration ----------
def cvz(a, n):
    # sum_{k>=0} (-1)^k a(k)
    d = (3 + mp.sqrt(8)) ** n
    d = (d + 1 / d) / 2
    b, c, s = mp.mpf(-1), -d, mp.mpf(0)
    for k in range(n):
        c = b - c
        s = s + c * a(k)
        b = (k + n) * (k - n) * b / ((k + mp.mpf(1) / 2) * (k + 1))
    return s / d


mp.mp.dps = 110
H = [mp.mpf(0)]
for k in range(1, 400):
    H.append(H[-1] + mp.mpf(1) / k)


def alt_H(p, n=160):
    # sum_{m>=1} (-1)^{m+1} H_m / m^p
    return cvz(lambda k: H[k + 1] / mp.mpf(k + 1) ** p, n)


def eta(m):
    return (1 - mp.mpf(2) ** (1 - m)) * mp.zeta(m)


z = mp.zeta
L2 = mp.log(2)
S3 = alt_H(3)
S3b = alt_H(3, 140)
classical = mp.mpf(11) / 4 * z(4) - mp.mpf(7) / 4 * z(3) * L2 + z(2) * L2 ** 2 / 2 - L2 ** 4 / 12 - 2 * mp.polylog(4, mp.mpf(1) / 2)
P_('sum (-1)^{n+1} H_n/n^3 (CVZ n=160) =', mp.nstr(S3, 50))
P_('  CVZ n=160 vs n=140 difference:', mp.nstr(S3 - S3b, 3))
P_('  minus classical closed form:', mp.nstr(S3 - classical, 3))
mp.mp.dps = 60
nsum_val = mp.nsum(lambda m: (-1) ** (int(m) + 1) * (mp.euler + mp.digamma(m + 1)) / m ** 3, [1, mp.inf])
P_('  nsum (dps 60) cross-check minus CVZ:', mp.nstr(nsum_val - S3, 3))
mp.mp.dps = 110

gamma12_claim = (-mp.mpf(15) / 32 * z(4) + mp.mpf(7) / 16 * z(3) * L2 - z(2) * L2 ** 2 / 8 + L2 ** 4 / 48
                 + mp.polylog(4, mp.mpf(1) / 2) / 2)

# gamma_{a,2}: Euler sum vs quadrature
mp.mp.dps = 60
def gq(a, b):
    f = lambda s: s ** a * mp.log1p(mp.exp(-2 * s)) ** b
    return mp.quad(f, [0, 1, 4, 12, 40, mp.inf]) / (mp.factorial(a) * mp.factorial(b))


G2 = {}
for a in range(0, 7):
    mp.mp.dps = 110
    es = mp.mpf(2) ** (-(a + 1)) * (eta(a + 3) - alt_H(a + 2))
    G2[a] = es
    mp.mp.dps = 60
    P_('gamma_{%d,2}: Euler sum - quadrature = %s' % (a, mp.nstr(es - gq(a, 2), 3)))
mp.mp.dps = 110
P_('gamma_12 (Euler sum) - claimed Li4 closed form =', mp.nstr(G2[1] - gamma12_claim, 3))

# PSLQ: odd a against pure zeta bases (no relation expected), even a (relation expected)
mp.mp.dps = 100
bases = {1: [z(4)], 2: [z(5), z(2) * z(3)], 3: [z(6), z(3) ** 2], 4: [z(7), z(2) * z(5), z(3) * z(4)],
         5: [z(8), z(3) * z(5), z(2) * z(3) ** 2]}
for a, B in bases.items():
    rel = mp.pslq([G2[a]] + B, maxcoeff=10 ** 12, maxsteps=10 ** 6)
    P_('PSLQ gamma_{%d,2} vs pure basis:' % a, rel)
rel = mp.pslq([G2[1], z(4), z(3) * L2, z(2) * L2 ** 2, L2 ** 4, mp.polylog(4, mp.mpf(1) / 2)], maxcoeff=10 ** 8,
              maxsteps=10 ** 6)
P_('PSLQ gamma_12 vs alternating weight-4 basis:', rel)

# ---------- direct F, Taylor coefficients by Cauchy contour ----------
# F_k(tau) := int_0^inf e^v v^k [(1+e^{-2v/tau})^tau - 1] dv = tau^{k+1} int_0^inf e^{tau s} s^k B(s,tau) ds
# (s-form = analytic continuation, analytic for Re tau < 2).
mp.mp.dps = 50
KS = [0, 1, 3, 5]


def Fk_all(tau):
    res = {}
    for k in KS:
        f = lambda s: mp.exp(tau * s) * s ** k * (mp.exp(tau * mp.log1p(mp.exp(-2 * s))) - 1)
        res[k] = tau ** (k + 1) * mp.quad(f, [0, 1, 3, 8, 20, 50, mp.inf])
    return res


# sanity: real tau, compare s-form with v-form for phi = e^v(1+v^3), tau = 0.3
tau0 = mp.mpf('0.3')
vform = mp.quad(lambda v: mp.exp(v) * (1 + v ** 3) * ((1 + mp.exp(-2 * v / tau0)) ** tau0 - 1),
                [0, tau0, 4 * tau0, 3, 10, 40, mp.inf])
sv = Fk_all(tau0)
P_('sanity: v-form F(0.3) for phi=e^v(1+v^3) minus s-form =', mp.nstr(vform - sv[0] - sv[3], 3))

N, r = 96, mp.mpf('0.75')
vals = {}
for m in range(N // 2 + 1):
    t = r * mp.expjpi(mp.mpf(2 * m) / N)
    vals[m] = Fk_all(t)
for m in range(N // 2 + 1, N):
    vals[m] = {k: mp.conj(vals[N - m][k]) for k in KS}


def coeff(k, n):
    s = mp.fsum(vals[m][k] * mp.expjpi(-mp.mpf(2 * m * n) / N) for m in range(N))
    return (s / N / r ** n).real


# Gamma closed form coefficients (psi = 1): 1 - Gamma(1-t/2)^2/Gamma(1-t)
gc = mp.taylor(lambda t: 1 - mp.gamma(1 - t / 2) ** 2 / mp.gamma(1 - t), 0, 10)
for n in range(0, 10):
    P_('psi=1: [tau^%d] contour - Gamma closed form = %s' % (n, mp.nstr(coeff(0, n) - gc[n], 3)))

mp.mp.dps = 50
for j in range(0, 3):
    k = 2 * j + 1
    n = 2 * j + 4
    # lower orders of F_k: zero below k+2, k! gamma_{k1} at k+2 (pure)
    low = [mp.nstr(coeff(k, m), 3) for m in range(0, k + 2)]
    first = coeff(k, k + 2) - mp.factorial(k) * eta(k + 2) / mp.mpf(2) ** (k + 1)
    pred_odd = mp.factorial(k + 1) * eta(k + 3) / mp.mpf(2) ** (k + 2) + mp.factorial(k) * G2[k]
    for omega in (mp.mpf(1), mp.mpf(3) / 7):
        Fn = coeff(0, n) + omega * coeff(k, n)
        pred = gc[n] + omega * pred_odd
        P_('j=%d omega=%s: F_%d contour = %s ; Prop E prediction = %s ; diff = %s' % (
            j, mp.nstr(omega, 5), n, mp.nstr(Fn, 32), mp.nstr(pred, 32), mp.nstr(Fn - pred, 3)))
    P_('   F_%d coefficients below tau^%d: %s ; [tau^%d] - k!eta(k+2)/2^(k+1) = %s' % (
        k, k + 2, low, k + 2, mp.nstr(first, 3)))
    # what the formula WITHOUT the middle term would give
    P_('   omission of middle term would be off by', mp.nstr(mp.factorial(k + 1) * eta(k + 3) / mp.mpf(2) ** (k + 2), 10))

# test example of the text: psi=(1+v)e^{-v^2}: F4 = gamma12 + pi^4/640
mp.mp.dps = 110
P_('text test: F4 predicted gamma12 + pi^4/640 =', mp.nstr(G2[1] + mp.pi ** 4 / 640, 45))
open(__file__.replace('.py', '_output.txt'), 'w').write('\n'.join(out) + '\n')
