"""RULQ script 1: coefficients of L_q(p)/8 ~ 1 + sum_k N_k h^{k+1} (h = 1/p) for rulers q = Q(h)/h^d.

Families (argument FAM):
  lam:L      q = L p            (Q = L, d = 1)
  lin:L:C    q = L p + C        (Q = L + C h, d = 1)
  sq         q = p^2            (Q = 1, d = 2)
Exact formula (s-form, see result.md Sec. 2):
  L/8 = 2^{-h} + h kappa int_0^oo G(kappa s,h) [(1+e^{-s})^eta - 1] ds,  kappa = h^d/(Q(1-h)), eta = h^d/Q,
  G(y,h) = e^{-hy}(1+e^{-y})^{-1-h} = sum_j g_j(h) y^j,
  N_k = (-log2)^{k+1}/(k+1)! + sum_{j,i} J_{i,j}/i! [h^k] (kappa^{j+1} g_j(h) eta^i).
For d = 1 and linear Q it also computes
  T_k   = n = 1 Laplace term  = sum_m m! L^{-m-1} [h^{k-2-m}] (a_m(h)/Q(h)),  A = e^{(1-beta)v}(1+e^v)^{-1-h},
          beta(h) = (Q(h)(1-h) - L)/h,
  S_k   = leading singular part (l = 0, all orders in h) at v0 = +-i pi,
  F_k   = two-term law  2Re{ e^{(1-b0)y0} L^{-k} [ (k-2)! y0^{1-k} + L (k-3)! y0^{2-k} (psi(k-2)+gamma-Log y0 - b1 y0 - Q1/L) ] }.
Usage: ru1_family.py FAM JFILE K DPS   -> writes ru1_N_<fam>_K<K>.txt and prints tables."""
import sys
import mpmath as mp

FAM, JF, K, DPS = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
mp.mp.dps = DPS
n = K + 1
fac = [mp.factorial(i) for i in range(2 * n + 2)]

# ---------- univariate truncated series helpers (length n)
def smul(a, b, m=n):
    return [mp.fsum(a[i] * b[k - i] for i in range(k + 1) if i < len(a) and k - i < len(b)) for k in range(m)]
def sinv(a, m=n):
    r = [mp.mpf(0)] * m; r[0] = 1 / a[0]
    for k in range(1, m):
        r[k] = -mp.fsum(a[i] * r[k - i] for i in range(1, min(k, len(a) - 1) + 1)) / a[0]
    return r
def sexp(v, m=n):
    f = [mp.exp(v[0])] + [mp.mpc(0) if isinstance(v[0], mp.mpc) else mp.mpf(0)] * (m - 1)
    for k in range(1, m):
        f[k] = mp.fsum(j * v[j] * f[k - j] for j in range(1, k + 1) if j < len(v)) / k
    return f
def slog(a, m=n):
    y = [mp.log(a[0])] + [mp.mpf(0)] * (m - 1)
    for k in range(1, m):
        y[k] = (k * a[k] - mp.fsum(j * y[j] * a[k - j] for j in range(1, k))) / (k * a[0])
    return y

# ---------- family
if FAM.startswith('lam:'):
    L = mp.mpf(mp.fraction(*[int(x) for x in FAM[4:].split('/')]) if '/' in FAM else FAM[4:]); C = mp.mpf(0); d = 1
elif FAM.startswith('lin:'):
    a, b = FAM[4:].split(':'); L = mp.mpf(a); C = mp.mpf(b); d = 1
elif FAM == 'sq':
    L = mp.mpf(1); C = mp.mpf(0); d = 2
else:
    raise SystemExit('unknown family')
Qs = [L, C] + [mp.mpf(0)] * (n - 2)
invQ = sinv(Qs)
hd = [mp.mpf(0)] * n; hd[d] = mp.mpf(1)
eta = smul(hd, invQ)
kappa = smul(eta, [mp.mpf(1)] * n)          # 1/(1-h) = sum h^i

# ---------- J table
J = {}
for line in open(JF):
    i, j, v = line.split(); J[(int(i), int(j))] = mp.mpf(v)

# ---------- g_j(h) = [y^j] e^{-hy}(1+e^{-y})^{-1-h}  (as in PIPDIV pd1)
lam2 = mp.log(2)
ey = [mp.mpf(1) / fac[i] for i in range(n + 1)]
onepe = [ey[0] + 1] + ey[1:]
inv = sinv(onepe, n + 1)
lser = [lam2] + [-inv[m - 1] / m for m in range(1, n)]     # l(y) = log(1+e^{-y})
X0 = [-c for c in lser]
X1 = [-lser[0], -lser[1] - 1] + [-c for c in lser[2:]]
E0 = sexp(X0)
G = [[mp.mpf(0)] * n for _ in range(n)]                     # G[j][r] = [y^j h^r]
P = E0[:]
for r in range(n):
    for j in range(n): G[j][r] = P[j] / fac[r]
    P = smul(P, X1)

# ---------- assemble N_k
etap = [None, eta]
for i in range(2, n): etap.append(smul(etap[-1], eta))
contrib = [mp.mpf(0)] * n
kp = kappa[:]                                               # kappa^{j+1}
for j in range(0, K):
    Pj = smul(kp, [G[j][r] for r in range(n)])
    W = [mp.mpf(0)] * n
    for i in range(1, K - j):
        if (i, j) not in J: continue
        c = J[(i, j)] / fac[i]
        e = etap[i]
        for t in range(n): W[t] += c * e[t]
    prod = smul(Pj, W)
    for k in range(n): contrib[k] += prod[k]
    kp = smul(kp, kappa)
N = [(-lam2) ** (k + 1) / fac[k + 1] + contrib[k] for k in range(n)]
tag = FAM.replace(':', '_').replace('/', 'o')
with open('ru1_N_%s_K%d.txt' % (tag, K), 'w') as f:
    for v in N: f.write(mp.nstr(v, DPS) + '\n')
print('family %s (lambda=%s, c=%s, d=%d), K=%d, dps=%d' % (FAM, L, C, d, K, DPS))
for k in range(6): print('N_%d = %s' % (k, mp.nstr(N[k], 30)))
if d != 1:
    print('k  N_k   |N_k|^(2/k) / (k/(2 e pi))')
    for k in range(4, K + 1):
        if N[k] != 0:
            print('%3d % .12e  %.6f' % (k, float(N[k]), float(abs(N[k]) ** (mp.mpf(2) / k) / (k / (2 * mp.e * mp.pi)))))
    raise SystemExit

# ---------- Laplace-form n = 1 term T_k
b0 = C - L; b1 = -C                                         # beta(h) = C - L - C h for linear Q
y0 = mp.mpc(0, mp.pi)
lp = slog([ey[0] + 1] + ey[1:])                             # log(1+e^v)
base = sexp([(1 - b0 if i == 1 else 0) - lp[i] for i in range(n)])   # e^{(1-b0)v}/(1+e^v)
X1v = [(-b1 if i == 1 else 0) - lp[i] for i in range(n)]    # h-coefficient of exponent: -b1 v - log(1+e^v)
Ar = [[mp.mpf(0)] * n for _ in range(n)]
P = base[:]
for r in range(n):
    for m in range(n): Ar[m][r] = P[m] / fac[r]
    P = smul(P, X1v)
T = [None] * n
for k in range(2, n):
    s = mp.mpf(0)
    for m in range(0, k - 1):
        am = smul([Ar[m][r] for r in range(n)], invQ, k - 1 - m)
        s += fac[m] * L ** (-m - 1) * am[k - 2 - m]
    T[k] = s
# ---------- leading singular part S_k (l = 0, exact in h)
w = C * y0 - mp.log(y0)
eh = sexp([mp.mpc(0), w] + [mp.mpc(0)] * (n - 2))
pre = mp.exp((1 - b0) * y0)
Sk = [None] * n
poch = [mp.mpf(1)] + [mp.mpf(0)] * (n - 1)                  # (1+h)_m as polynomial
pochs = [poch[:]]
for m in range(1, n):
    poch = smul(poch, [mp.mpf(m), mp.mpf(1)]); pochs.append(poch[:])
ehQ = smul(eh, invQ)
for k in range(2, n):
    s = mp.mpc(0)
    for m in range(0, k - 1):
        ser = smul(ehQ, pochs[m], k - 1 - m)
        s += L ** (-m - 1) * y0 ** (-1 - m) * ser[k - 2 - m]
    Sk[k] = 2 * mp.re(pre * s)
# ---------- two-term law
def F(k):
    br = mp.digamma(k - 2) + mp.euler - mp.log(y0) - b1 * y0 - C / L
    z = fac[k - 2] * y0 ** (1 - k) + L * fac[k - 3] * y0 ** (2 - k) * br
    return 2 * mp.re(pre * L ** (-k) * z)
print('phase: leading N_k ~ -2(k-2)!/(pi^{k-1} lambda^k) cos(pi(lambda-c) - pi(k-1)/2)')
print(' k  N_k                 nu_k=N_k/scale  F_k/scale     T_k/N_k-1     (N-S)/scale*k^2   (N-F)/scale*k^2/log^2k')
for k in list(range(4, 21)) + list(range(21, K + 1, 4)) + [K - 1, K]:
    if k < 4 or k > K: continue
    scale = 2 * fac[k - 2] / (mp.pi ** (k - 1) * L ** k)
    print('%3d % .12e % .9f % .9f % .3e % .4f % .4f' % (k, float(N[k]), float(N[k] / scale), float(F(k) / scale),
          float(T[k] / N[k] - 1) if N[k] != 0 else 0.0, float((N[k] - Sk[k]) / scale * k**2),
          float((N[k] - F(k)) / scale * k**2 / mp.log(k)**2)))
