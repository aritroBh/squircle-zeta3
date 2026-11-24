import json
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, nstr, findroot
import numpy as np
mp.dps = 200
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
def peval(c, t):
    v = c[0]
    for k in range(1,len(c)): v += c[k]*t**k
    return v
num90, den90 = pade_coeffs(b1[:181], 90, 90)
c = [complex(den90[k]) for k in range(len(den90)-1,-1,-1)]
rts = np.roots(c)
f = lambda t: peval(den90, t)
ok, fail = [], []
for r in rts:
    try:
        rt = findroot(f, mpc(r.real, r.imag), solver='mnewton', tol=mpf('1e-25'), maxsteps=60, verify=True)
        ok.append(rt)
    except Exception:
        fail.append(r)
print("refined: %d, failed: %d" % (len(ok), len(fail)))
prs = sorted(ok, key=lambda r: abs(r))
print("\nPB1[90/90] TRUE (refined) poles nearest origin:")
for r in prs[:10]:
    print("   t = %s   |t| = %s   |den| = %s" % (nstr(r,18), nstr(abs(r),14), nstr(abs(peval(den90,r)),2)))
