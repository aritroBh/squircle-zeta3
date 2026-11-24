import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, nstr
from math import comb
t0=time.time()
mp.dps = 250
f450 = [mpf(v) for v in json.load(open('step2_num450_N180.json'))['f450']]
sq2 = sqrt(2)
f1 = [f450[n] - sq2*((-1)**n) for n in range(len(f450))]
b1 = [f1[n]/factorial(n) for n in range(len(f1))]
NMAX = 180
ipi = mpc(0,1)*pi
# d_k = sum_{m=0}^{floor(k/2)} b1_{k-2m} (i pi)^{k-2m} (-1)^m C(k-m-1, m)   (k>=1), d_0 = b1_0
d = [mpc(0)]*(NMAX+1)
d[0] = b1[0]
for k in range(1, NMAX+1):
    s = mpc(0)
    for m in range(0, k//2+1):
        n = k-2*m
        if n == 0:
            continue  # t^0 contributes only to k=0
        s += b1[n] * ipi**n * (-1)**m * comb(k-m-1, m)
    d[k] = s
print("built d-series %.1fs; sample:"%(time.time()-t0))
for k in [1,2,10,50,100,150,180]:
    print("  d_%d = %s   |d_k|=%s" % (k, nstr(d[k],10), nstr(abs(d[k]),4)))
# radius check
print("|d_k|^(1/k):", [nstr(abs(d[k])**(1/mpf(k)),5) for k in [90,120,150,180]])
def pade_coeffs(a, L, M):
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            k2 = L + i - j
            A[i, j] = a[k2] if k2 >= 0 else mpf(0)
        rhs[i] = -a[L + 1 + i]
    q = (A**-1)*rhs
    den = [mpf(1)] + [q[j] for j in range(M)]
    num = [sum(den[j]*a[n-j] for j in range(0, min(M, n)+1)) for n in range(L+1)]
    return num, den
numW85, denW85 = pade_coeffs(d[:171], 85, 85)
numW90, denW90 = pade_coeffs(d[:181], 90, 90)
json.dump({'num':[nstr(v,200) for v in numW90],'den':[nstr(v,200) for v in denW90],
           'num85':[nstr(v,200) for v in numW85],'den85':[nstr(v,200) for v in denW85],
           'd':[nstr(v,200) for v in d]}, open('step_12_wpade.json','w'))
print("w-Pade built+saved %.1fs"%(time.time()-t0))
# poles of w-pade (numpy, but coefficient range is mild here)
import numpy as np
for tag,(num,den) in [('85/85',(numW85,denW85)),('90/90',(numW90,denW90))]:
    c = [complex(den[k]) for k in range(len(den)-1,-1,-1)]
    rts = np.roots(c)
    near1 = sorted(rts, key=lambda r: abs(r-1))[:6]
    print(f"w-Pade[{tag}] poles nearest w=1:")
    for r in near1:
        print("     %+.8f %+.8fi   |w|=%.8f  |w-1|=%.6f" % (r.real, r.imag, abs(r), abs(r-1)))
