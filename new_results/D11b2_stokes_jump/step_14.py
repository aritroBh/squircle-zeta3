import time
from mpmath import mp, mpf, mpc, pi, sqrt, nstr
import numpy as np
mp.dps = 200
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
def pat(n): return [-1,1,1,-1][n%4]
tp = 2/pi
# MODEL 2-cut: a_n = pat(n) n^{-3/2} (2/pi)^n  (exact D11 leading asymptotic, both tips, correct phases)
a2 = [mpf(0)] + [pat(n)*mpf(n)**mpf('-1.5')*tp**n for n in range(1,181)]
# MODEL 1-cut: a_n = i^n n^{-3/2} (2/pi)^n -> Li_{3/2}(2it/pi): single cut at -i pi/2 ray
a1 = [mpf(0)] + [mpc(0,1)**n * mpf(n)**mpf('-1.5')*tp**n for n in range(1,181)]
for tag, a in [('2-cut model (pat phases)', a2), ('1-cut model (Li_{3/2}(2it/pi))', a1)]:
    for N in [70, 90]:
        num, den = pade_coeffs(a[:2*N+1], N, N)
        c = [complex(den[k]) for k in range(len(den)-1,-1,-1)]
        rts = np.roots(c)
        prs = sorted(rts, key=lambda r: abs(r))[:4]
        print(f"{tag}  [{N}/{N}] innermost poles:")
        for r in prs:
            print("     %+.10f %+.10fi   |t|=%.10f" % (r.real, r.imag, abs(r)))
    print()
