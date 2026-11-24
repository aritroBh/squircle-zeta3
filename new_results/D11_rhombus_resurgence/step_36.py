# Pole structure via mpmath polyroots (high precision) + proper Pade series-residual check
from mpmath import mp, mpf, mpc, pi, factorial, polyroots
import json
mp.dps = 100
Fcoef = [mpf(v) for v in json.load(open('step2_num200.json'))['f200']][:62]
b = [Fcoef[n]/factorial(n) for n in range(61)]
def pade_coeffs(a, L, M):
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            A[i, j] = a[L - M + i + j + 1] if L - M + i + j + 1 >= 0 else mpf(0)
        rhs[i] = -a[L + i + 1]
    q = (A**-1)*rhs
    den = [mpf(1)] + [q[i] for i in range(M)]
    num = [mpf(0)]*(L+1)
    for n in range(L+1):
        num[n] = sum(den[j]*a[n-j] for j in range(0, min(M, n)+1))
    return num, den
num, den = pade_coeffs(b, 30, 30)
# residual: series of num/den vs b
def series_div(num, den, N):
    out = [mpf(0)]*(N+1)
    for n in range(N+1):
        s = num[n] if n < len(num) else mpf(0)
        out[n] = (s - sum(den[j]*out[n-j] for j in range(1, min(len(den)-1, n)+1)))/den[0]
    return out
rec = series_div(num, den, 60)
print("max |reconstructed b_n - b_n|, n=0..60:", mp.nstr(max(abs(rec[n]-b[n]) for n in range(61)), 3))
print("max |b_n| ~", mp.nstr(max(abs(v) for v in b), 3))
roots = polyroots([den[k] for k in range(30, -1, -1)], maxsteps=200)
roots_sorted = sorted(roots, key=lambda r: abs(r))
print("\npoles of Pade[30/30] by modulus (first 8):")
for r in roots_sorted[:8]:
    print("   ", mp.nstr(r, 20), "   |t| =", mp.nstr(abs(r), 10))
