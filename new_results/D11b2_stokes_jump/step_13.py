import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, nstr
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

# 1) poles of the SAVED pade1_9090.json (built at dps=400 in D11b)
pj = json.load(open('pade1_9090.json'))
numS = [mpf(v) for v in pj['num']]; denS = [mpf(v) for v in pj['den']]
cS = [complex(denS[k]) for k in range(len(denS)-1,-1,-1)]
rtsS = np.roots(cS)
prsS = sorted(rtsS, key=lambda r: abs(r))[:6]
print("SAVED pade1_9090.json (dps400 build) innermost poles:")
for r in prsS:
    print("   %+.10f %+.10fi   |t|=%.10f" % (r.real, r.imag, abs(r)))

# 2) rebuild [90/90] at dps=350 and compare
mp.dps = 350
t0=time.time()
num90hi, den90hi = pade_coeffs(b1[:181], 90, 90)
print("rebuilt [90/90] at dps=350 (%.1fs)"%(time.time()-t0))
chi = [complex(den90hi[k]) for k in range(len(den90hi)-1,-1,-1)]
rtshi = np.roots(chi)
prshi = sorted(rtshi, key=lambda r: abs(r))[:6]
print("dps350-build [90/90] innermost poles:")
for r in prshi:
    print("   %+.10f %+.10fi   |t|=%.10f" % (r.real, r.imag, abs(r)))

# difference between the two builds' denominators (coefficient scale)
mp.dps = 120
dmax = max(abs(den90hi[k]-denS[k]) for k in range(91))
smax = max(abs(denS[k]) for k in range(91))
print("\nmax |den diff| between dps350 and saved builds = %s ; max|den| = %s" % (nstr(dmax,3), nstr(smax,3)))
