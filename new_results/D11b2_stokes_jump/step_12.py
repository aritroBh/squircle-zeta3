import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, nstr
import numpy as np
mp.dps = 120
f450 = [mpf(v) for v in json.load(open('step2_num450_N180.json'))['f450']]
sq2 = sqrt(2)
f1 = [f450[n] - sq2*((-1)**n) for n in range(len(f450))]
b1 = [f1[n]/factorial(n) for n in range(len(f1))]
b  = [f450[n]/factorial(n) for n in range(len(f450))]
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

mp.dps = 200
num90, den90 = pade_coeffs(b1[:181], 90, 90)
# verify numpy roots are true roots: evaluate den90 at them with mp
import numpy.polynomial.polynomial as P
c = [complex(den90[k]) for k in range(len(den90)-1,-1,-1)]
rts = np.roots(c)
prs = sorted(rts, key=lambda r: abs(r))[:4]
print("verify den90 at numpy roots (should be ~0):")
for r in prs:
    val = peval(den90, mpc(r.real, r.imag))
    # residual of the actual rational function near there: residue estimate
    print("  root %+.12f%+.12fi  den=%s" % (r.real, r.imag, nstr(abs(val),3)))

# residues at the inner poles: Res = num/den' 
def pder(c):
    return [c[k]*k for k in range(1,len(c))]
print("\nresidues num/den' at inner poles:")
dden = pder(den90)
for r in prs:
    t = mpc(r.real, r.imag)
    res = peval(num90, t)/peval(dden, t)
    print("  pole %+.10f%+.10fi  |Res| = %s" % (r.real, r.imag, nstr(abs(res),6)))

# ---- Padé of the UNSUBTRACTED B: does it also show inner poles? ----
numB, denB = pade_coeffs(b[:181], 90, 90)
cB = [complex(denB[k]) for k in range(len(denB)-1,-1,-1)]
rtsB = np.roots(cB)
prsB = sorted(rtsB, key=lambda r: abs(r))[:6]
print("\nPB (unsubtracted) [90/90] poles nearest origin:")
for r in prsB:
    print("   %+.10f %+.10fi   |t|=%.10f" % (r.real, r.imag, abs(r)))
