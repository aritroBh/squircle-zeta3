"""Referee script 5: N_k for k <= KMAX by numerical power-series arithmetic in h at each
quadrature node s (Gauss-Legendre on unit intervals), to test the growth claim
|N_k| ~ 2 (k-2)!/pi^(k-1).  Independent of the sympy route of r1 (checked against it for k<=8)."""
import sys, pickle
import mpmath as mp
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
mp.mp.dps = 70
K = KMAX + 2
lam = mp.log(2)

def sexp0(v):  # exp of series with v[0]=0
    f = [mp.mpf(0)] * (K + 1); f[0] = mp.mpf(1)
    for n in range(1, K + 1):
        f[n] = mp.fsum(k * v[k] * f[n - k] for k in range(1, n + 1)) / n
    return f
def slog1p(q):  # log(1+q), q[0]=0
    y = [mp.mpf(0)] * (K + 1)
    for n in range(1, K + 1):
        y[n] = (n * q[n] - mp.fsum(k * y[k] * q[n - k] for k in range(1, n))) / n
    return y
def smul(a, b):
    return [mp.fsum(a[i] * b[n - i] for i in range(n + 1)) for n in range(K + 1)]

A = [mp.mpf(0)] + [mp.mpf(1)] * K            # h/(1-h)
hA = [mp.mpf(0)] + A[:K]                     # h^2/(1-h)
def coeffs(s):
    u = [s * x for x in A]
    em = sexp0([-x for x in u])
    q = [mp.mpf(0)] + [x / 2 for x in em[1:]]
    g = slog1p(q)                            # log(1+e^{-u}) - lam
    onehg = [g[n] + (g[n - 1] if n > 0 else 0) for n in range(K + 1)]   # (1+h) g
    E1 = [-s * hA[n] - onehg[n] - (lam if n == 1 else 0) for n in range(K + 1)]
    w = [x / 2 for x in sexp0(E1)]
    L = mp.log1p(mp.exp(-s))
    br = [mp.mpf(0)] + [L**j / mp.factorial(j) for j in range(1, K + 1)]
    return smul(hA, smul(w, br))

xs, ws = zip(*[(x, w) for x, w in zip(*mp.gauss_quadrature(40, 'legendre'))]) if hasattr(mp, 'gauss_quadrature') else (None, None)
if xs is None:
    import mpmath.calculus.quadrature as Q
    gl = Q.GaussLegendre(mp)
    nodes = gl.calc_nodes(6, mp.prec)  # degree 6 -> 3*2^5 = 96 nodes on [-1,1]
    xs, ws = zip(*nodes)
tot = [mp.mpf(0)] * (K + 1)
SMAX = 3 * KMAX + 160
for a0 in range(0, SMAX, 2):
    for x, w in zip(xs, ws):
        s = a0 + 1 + x
        c = coeffs(s)
        for n in range(K + 1): tot[n] += w * c[n]
N = [(-lam)**(k + 1) / mp.factorial(k + 1) + tot[k + 1] for k in range(KMAX + 1)]
ref = pickle.load(open('r1_Nnum.pkl', 'rb'))
for k in range(min(9, KMAX + 1)):
    print('N_%d check vs r1: %s' % (k, mp.nstr(N[k] - mp.mpf(ref[k]), 3)))
for k in range(2, KMAX + 1):
    norm = N[k] * mp.pi**(k - 1) / (2 * mp.factorial(k - 2))
    print('k=%2d  N_k=% .12e   N_k pi^(k-1)/(2(k-2)!) = % .6f   |N_k|^(1/k)=%.4f  k/(e pi)=%.4f' % (
        k, float(N[k]), float(norm), float(abs(N[k])**(mp.mpf(1) / k)), k / float(mp.e * mp.pi)))
