from mpmath import mp, mpf, mpc, pi, factorial, polyroots
import json
mp.dps = 110
Fcoef = [mpf(v) for v in json.load(open('step2_num200.json'))['f200']][:62]
b = [Fcoef[n]/factorial(n) for n in range(61)]
def pade_coeffs(a, L, M):
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            k = L + i - j
            A[i, j] = a[k] if k >= 0 else mpf(0)
        rhs[i] = -a[L + 1 + i]
    q = (A**-1)*rhs
    den = [mpf(1)] + [q[j] for j in range(M)]
    num = [sum(den[j]*a[n-j] for j in range(0, min(M, n)+1)) for n in range(L+1)]
    return num, den
num, den = pade_coeffs(b, 30, 30)
def series_div(num, den, N):
    out = [mpf(0)]*(N+1)
    for n in range(N+1):
        s = num[n] if n < len(num) else mpf(0)
        out[n] = (s - sum(den[j]*out[n-j] for j in range(1, min(len(den)-1, n)+1)))/den[0]
    return out
rec = series_div(num, den, 60)
print("max |reconstructed b_n - b_n| n<=60:", mp.nstr(max(abs(rec[n]-b[n]) for n in range(61)), 3))
roots = polyroots([den[k] for k in range(30, -1, -1)], maxsteps=300)
rs = sorted(roots, key=lambda r: abs(r))
print("nearest poles:")
for r in rs[:6]:
    print("   ", mp.nstr(r, 18), "  |t| =", mp.nstr(abs(r), 12))
json.dump({'num': [mp.nstr(v, 95) for v in num], 'den': [mp.nstr(v, 95) for v in den]}, open('step3_pade.json','w'))
