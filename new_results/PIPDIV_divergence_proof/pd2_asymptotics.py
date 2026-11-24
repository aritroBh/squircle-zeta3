"""PIPDIV script 2: test the proved structure of N_k against the numbers of pd1.

  N_k = (-log2)^{k+1}/(k+1)! + T_k + U_k,
  T_k = sum_{m=0}^{k-2} m! [y^m h^{k-2-m}] A(y,h),  A = e^{2y} (1+e^y)^{-1-h}   (the n = 1 term)
  U_k = rest (n >= 2 terms; proved O(2^k + k^2 (k-2)! 6^{-k}))
Singular expansion of A at y0 = +-i pi:  A = (y0-y)^{-1-h} Phi(y-y0,h),  Phi = e^{2u} ((e^u-1)/u)^{-1-h},
  m! [y^m] u^l (y0-y)^{-1-h} = (-1)^l y0^{l-1-h-m} (1+h-l)_m.
T^{(L)}_k = sum over y0, l < L of sum_m [h^{k-2-m}] phi_l(h) (-1)^l y0^{l-1-h-m} (1+h-l)_m.
Two-term law F_k = 2 Re sum_{y0=i pi} [ (k-2)! y0^{1-k} + (k-3)! y0^{2-k} (psi(k-2) + gamma - Log y0) ].
Usage: pd2_asymptotics.py NFILE K"""
import sys
import mpmath as mp
mp.mp.dps = 60
NF = sys.argv[1]; K = int(sys.argv[2])
N = [mp.mpf(l.strip()) for l in open(NF)][:K + 1]
n = K + 2
lam = mp.log(2)

def smul(a, b):
    m = min(len(a), len(b))
    return [mp.fsum(a[i] * b[j - i] for i in range(j + 1)) for j in range(m)]
def sexp(v):
    f = [mp.exp(v[0])] + [mp.mpf(0)] * (len(v) - 1)
    for m in range(1, len(v)):
        f[m] = mp.fsum(k * v[k] * f[m - k] for k in range(1, m + 1)) / m
    return f
def slog(a):  # log of series with a[0] > 0
    y = [mp.log(a[0])] + [mp.mpf(0)] * (len(a) - 1)
    for m in range(1, len(a)):
        y[m] = (m * a[m] - mp.fsum(k * y[k] * a[m - k] for k in range(1, m))) / (m * a[0])
    return y
fac = [mp.factorial(i) for i in range(n + 2)]

# --- A(y,h): Ar[m][r]
ey = [1 / fac[i] for i in range(n)]
lt = slog([ey[0] + 1] + ey[1:])            # log(1+e^y)
base = sexp([(2 if i == 1 else 0) - lt[i] for i in range(n)])   # e^{2y}/(1+e^y)
Ar = [[mp.mpf(0)] * n for _ in range(n)]
P = base[:]
mlt = [-c for c in lt]
for r in range(n):
    for m in range(n): Ar[m][r] = P[m] / fac[r]
    P = smul(P, mlt)
T = [mp.fsum(fac[m] * Ar[m][k - 2 - m] for m in range(k - 1)) if k >= 2 else None for k in range(K + 1)]

# --- Phi(u,h) coefficients phi[l][r]
LMAX = 5
phis = slog([1 / fac[i + 1] for i in range(n)])  # log((e^u-1)/u)
Phi = [[mp.mpf(0)] * n for _ in range(LMAX)]
B0 = sexp([(2 if i == 1 else 0) - phis[i] for i in range(n)])  # e^{2u} phi^{-1}
P = B0[:]
mph = [-c for c in phis]
for r in range(n):
    for l in range(LMAX): Phi[l][r] = P[l] / fac[r]
    P = smul(P, mph)

def poch_series(a0, m, deg):  # (a0 + h)_m as h-series, truncated at deg
    p = [mp.mpf(1)] + [mp.mpf(0)] * deg
    for i in range(m):
        c = a0 + i
        p = [c * p[0]] + [c * p[j] + p[j - 1] for j in range(1, deg + 1)]
    return p
def TL(k, L):
    tot = mp.mpc(0)
    for y0 in (mp.mpc(0, mp.pi), mp.mpc(0, -mp.pi)):
        ly = mp.log(y0)
        for l in range(L):
            for m in range(k - 1):
                d = k - 2 - m
                ex = [(-ly) ** i / fac[i] for i in range(d + 1)]
                s = smul(smul(ex, Phi[l][:d + 1]), poch_series(1 - l, m, d))
                tot += (-1) ** l * y0 ** (l - 1 - m) * s[d]
    return tot.real
def F2(k):
    y0 = mp.mpc(0, mp.pi)
    v = fac[k - 2] * y0 ** (1 - k) + fac[k - 3] * y0 ** (2 - k) * (mp.psi(0, k - 2) + mp.euler - mp.log(y0))
    return 2 * v.real

print('phi_1(h) = %s + %s h  (expect 3/2 - h/2)' % (mp.nstr(Phi[1][0], 15), mp.nstr(Phi[1][1], 15)))
print('\nTable 1: U_k = N_k - (-log2)^{k+1}/(k+1)! - T_k  versus the scale (k-2)!/(2pi)^{k-1}')
for k in range(2, K + 1):
    U = N[k] - (-lam) ** (k + 1) / fac[k + 1] - T[k]
    print('k=%2d  T_k/N_k=% .12f  U_k (2pi)^{k-1}/(k-2)! = % .5f' % (k, float(T[k] / N[k]), float(U * (2 * mp.pi) ** (k - 1) / fac[k - 2])))

def F3(k):  # three-term law: adds the i=2 term of P_k and the i=1 term of Q_k
    y0 = mp.mpc(0, mp.pi); Ly = mp.log(y0); X = mp.psi(0, k - 3) + mp.euler
    v = fac[k - 4] * y0 ** (3 - k) * ((X - Ly) ** 2 / 2 + (mp.psi(1, k - 3) - mp.zeta(2)) / 2 - mp.mpf(3) / 2)
    return F2(k) + 2 * v.real

print('\nTable 2: normalised nu_k = N_k pi^{k-1}/(2(k-2)!), two- and three-term laws, scaled residuals,')
print('         and |T_k - T^(L)_k|/scale for L = 1..4 (scale = 2(k-2)!/pi^{k-1})')
print(' k   nu_k           2-term          3-term          r2=(nu-nu2)k^2/log^2k  r3=(nu-nu3)k^3/log^3k   L=1      L=2      L=3      L=4')
for k in range(5, K + 1):
    sc = 2 * fac[k - 2] / mp.pi ** (k - 1)
    nu = N[k] / sc
    f2 = F2(k) / sc; f3 = F3(k) / sc
    r2 = (nu - f2) * k ** 2 / mp.log(k) ** 2
    r3 = (nu - f3) * k ** 3 / mp.log(k) ** 3
    errs = [float(abs(T[k] - TL(k, L)) / sc) for L in (1, 2, 3, 4)]
    print('%3d % .10f % .10f % .10f  % .4f  % .4f   %.2e %.2e %.2e %.2e' % (k, float(nu), float(f2), float(f3), float(r2), float(r3), *errs))
