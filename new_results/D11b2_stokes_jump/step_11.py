import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, log, nstr
import numpy as np
mp.dps = 120
f450 = [mpf(v) for v in json.load(open('step2_num450_N180.json'))['f450']]
sq2 = sqrt(2)
f1 = [f450[n] - sq2*((-1)**n) for n in range(len(f450))]
b1 = [f1[n]/factorial(n) for n in range(len(f1))]
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
t0=time.time(); mp.dps=200
num70, den70 = pade_coeffs(b1[:141], 70, 70)
num90, den90 = pade_coeffs(b1[:181], 90, 90)
json.dump({'num':[nstr(v,180) for v in num90],'den':[nstr(v,180) for v in den90]}, open('step_09_pade90.json','w'))
json.dump({'num':[nstr(v,180) for v in num70],'den':[nstr(v,180) for v in den70]}, open('step_09_pade70.json','w'))
mp.dps=60
def roots_np(coef):
    c = [complex(coef[k]) for k in range(len(coef)-1,-1,-1)]
    return np.roots(c)
for tag,(num,den) in [('70/70',(num70,den70)),('90/90',(num90,den90))]:
    pr = roots_np(den); zr = roots_np(num)
    prs = sorted(pr, key=lambda r: abs(r))
    print(f"==== PB1[{tag}] poles nearest origin ====")
    for r in prs[:8]:
        dz = min(abs(r-z) for z in zr)
        print("   pole %+.10f %+.10fi   |t|=%.10f   d_zero=%.2e" % (r.real, r.imag, abs(r), dz))
print("done %.1fs"%(time.time()-t0))
