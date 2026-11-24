import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, nstr, findroot
import numpy as np
t0=time.time()
mp.dps = 80
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
num70, den70 = pade_coeffs(b1[:141], 70, 70)
c = [complex(den70[k]) for k in range(len(den70)-1,-1,-1)]
rts = np.roots(c)
cand = [r for r in rts if abs(r) < 1.95]
res = []
for r in cand:
    try:
        rt = findroot(lambda t: peval(den70,t), mpc(r.real, r.imag), solver='mnewton', tol=mpf('1e-15'), maxsteps=20, verify=False)
        res.append(rt)
    except Exception: pass
res = sorted(set([mpc(r) for r in res]), key=lambda r: abs(r))
print("PB1[70/70] refined innermost TRUE poles (dps80):")
seen=[]
for r in res[:8]:
    print("   t=%s  |t|=%s  |den|=%s" % (nstr(r,15), nstr(abs(r),12), nstr(abs(peval(den70,r)),2)))
print("pi/2 = %s ; smallest |t| = %s" % (nstr(pi/2,15), nstr(min(abs(r) for r in res),12)))
print("(%.1fs)"%(time.time()-t0))
