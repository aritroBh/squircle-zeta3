import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, nstr, findroot
import numpy as np
mp.dps = 100
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
t0=time.time()
num90, den90 = pade_coeffs(b1[:181], 90, 90)
print("built %.1fs"%(time.time()-t0))
c = [complex(den90[k]) for k in range(len(den90)-1,-1,-1)]
rts = np.roots(c)
# only candidates with |t| < 2.5 (inner region) 
cand = [r for r in rts if abs(r) < 2.5]
print("numpy candidates inside |t|<2.5: %d of %d" % (len(cand), len(rts)))
ok = []
for r in cand:
    try:
        rt = findroot(lambda t: peval(den90,t), mpc(r.real, r.imag), solver='mnewton', tol=mpf('1e-20'), maxsteps=25, verify=False)
        moved = abs(rt - mpc(r.real,r.imag))
        ok.append((rt, moved, abs(peval(den90,rt))))
    except Exception as e:
        ok.append((None, None, None))
print("refined inner poles:")
for (rt, moved, resid) in sorted([x for x in ok if x[0] is not None], key=lambda x: abs(x[0])):
    print("   t=%s  |t|=%s  moved=%s  |den|=%s" % (nstr(rt,16), nstr(abs(rt),12), nstr(moved,2), nstr(resid,2)))
print("%.1fs"%(time.time()-t0))
