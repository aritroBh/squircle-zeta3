"""PIPDIV script 1: the coefficients N_k (pi_p/4 = 1 + sum_k N_k h^{k+1}, h = 1/p) for k <= K.

Route (independent of the referee's r5, which integrates h-series at quadrature nodes):
  N_k = (-log2)^{k+1}/(k+1)! + sum_{j=0}^{k-2} j! [h^k]( kappa^{j+1} g_j(h) B_j(h) ),
  kappa = h/(1-h),  g_j(h) = [y^j] e^{-hy} (1+e^{-y})^{-1-h},
  B_j(h) = sum_{r>=1} h^r J_{r,j} / (r! j!),  J_{r,j} = int_0^oo s^j log^r(1+e^{-s}) ds.
J table by composite Gauss-Legendre on [0, S] (unit-2 panels).  Usage: pd1_Nk.py K DPS NODEDEG S
Writes pd1_N_K{K}_dps{DPS}.txt (one N_k per line, full precision)."""
import sys, time
import mpmath as mp
import mpmath.calculus.quadrature as Q

K = int(sys.argv[1]); mp.mp.dps = int(sys.argv[2]); DEG = int(sys.argv[3]); S = int(sys.argv[4])
t0 = time.time()
lam = mp.log(2)

# ---- J table: J[r][j] for r>=1, r+j <= K-1
xs, ws = zip(*Q.GaussLegendre(mp.mp).calc_nodes(DEG, mp.mp.prec))
J = [[mp.mpf(0)] * K for _ in range(K + 1)]
for a0 in range(0, S, 2):
    for x, w in zip(xs, ws):
        s = a0 + 1 + x
        L = mp.log1p(mp.exp(-s))
        sp = [mp.mpf(1)]
        for j in range(1, K): sp.append(sp[-1] * s)
        Lr = w
        for r in range(1, K):
            Lr = Lr * L
            row = J[r]
            for j in range(0, K - r):
                row[j] += Lr * sp[j]
print('J table done %.0fs' % (time.time() - t0), flush=True)

# ---- y-series of G(y,h) = exp(X0 + h X1), X0 = -l(y), X1 = -l(y) - y, l = log(1+e^{-y})
def smul(a, b, n):
    return [mp.fsum(a[i] * b[m - i] for i in range(m + 1)) for m in range(n)]
def sexp(v, n):  # exp of series (any constant term)
    f = [mp.exp(v[0])] + [mp.mpf(0)] * (n - 1)
    for m in range(1, n):
        f[m] = mp.fsum(k * v[k] * f[m - k] for k in range(1, m + 1)) / m
    return f
n = K
# l(y) = log2 - y/2 + log cosh(y/2); Taylor of l via derivative: l' = -1/(1+e^{y})
ey = [mp.mpf(1) / mp.factorial(i) for i in range(n + 1)]
onepe = [ey[0] + 1] + ey[1:]
inv = [mp.mpf(0)] * (n + 1); inv[0] = 1 / onepe[0]
for m in range(1, n + 1):
    inv[m] = -mp.fsum(onepe[i] * inv[m - i] for i in range(1, m + 1)) / onepe[0]
lser = [lam] + [-inv[m - 1] / m for m in range(1, n)]
X0 = [-c for c in lser]
X1 = [-lser[0], -lser[1] - 1] + [-c for c in lser[2:]]
E0 = sexp(X0, n)
# G[j][r] = [y^j h^r]
G = [[mp.mpf(0)] * n for _ in range(n)]
P = E0[:]
for r in range(n):
    for j in range(n): G[j][r] = P[j] / mp.factorial(r)
    P = smul(P, X1, n)

# ---- assemble
def hser_mul(a, b, n):
    return [mp.fsum(a[i] * b[m - i] for i in range(m + 1)) for m in range(n)]
N = []
for k in range(K + 1):
    tot = (-lam) ** (k + 1) / mp.factorial(k + 1)
    for j in range(0, k - 1):
        # need [h^{k-j-2}] of (1-h)^{-j-1} g_j(h) b_j(h), b_j = B_j / h
        d = k - j - 2
        oneh = [mp.binomial(j + i, i) for i in range(d + 1)]
        g = [G[j][r] for r in range(d + 1)]
        b = [J[r + 1][j] / (mp.factorial(r + 1) * mp.factorial(j)) for r in range(d + 1)]
        prod = hser_mul(hser_mul(oneh, g, d + 1), b, d + 1)
        tot += mp.factorial(j) * prod[d]
    N.append(tot)
with open('pd1_N_K%d_dps%d_deg%d_S%d.txt' % (K, mp.mp.dps, DEG, S), 'w') as f:
    for v in N: f.write(mp.nstr(v, mp.mp.dps) + '\n')
z3 = mp.zeta(3)
cf = [-lam, lam**2 / 2, mp.pi**2 / 24 - lam**3 / 6,
      lam**4 / 24 + mp.pi**2 / 24 * (1 - lam) + z3 / 4,
      -lam**5 / 120 + mp.pi**2 / 24 * (1 - lam + lam**2 / 2) + z3 / 4 * (1 - lam) + mp.pi**4 / 2880,
      lam**6 / 720 + mp.pi**2 / 24 * (1 - lam + lam**2 / 2 - lam**3 / 6) + z3 / 4 * (1 - lam + lam**2 / 2)
      + mp.pi**4 / 2880 * (1 - lam) - 7 * mp.pi**4 / 1920 - z3 * mp.pi**2 / 48 + mp.zeta(5) / 8]
for k in range(6):
    print('N_%d - closed form = %s' % (k, mp.nstr(N[k] - cf[k], 3)))
for k in range(2, K + 1):
    print('k=%2d N_k=% .15e  norm=% .10f' % (k, float(N[k]), float(N[k] * mp.pi**(k - 1) / (2 * mp.factorial(k - 2)))))
print('total %.0fs' % (time.time() - t0))
