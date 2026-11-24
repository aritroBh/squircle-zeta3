import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, log10, nstr, cos
mp.dps = 120
t0=time.time()
f450 = [mpf(v) for v in json.load(open('step2_num450_N180.json'))['f450']]
print("loaded f_0..f_%d" % (len(f450)-1))
# subtract the geometric part sqrt2*(-1)^n  -> f1 (purely factorial-divergent part)
sq2 = sqrt(2)
f1 = [f450[n] - sq2*((-1)**n) for n in range(len(f450))]
b  = [f450[n]/factorial(n) for n in range(len(f450))]   # Borel coeffs of B
b1 = [f1[n]/factorial(n) for n in range(len(f1))]        # Borel coeffs of B1
def pat(n): return [-1,1,1,-1][n%4]

# ---- root test on b_n ----
print("\nroot test  |b_n|^(-1/n)  (should -> pi/2 = %s)" % nstr(pi/2,15))
for n in [20,40,60,80,100,120,140,160,180]:
    r = abs(b[n])**(-mpf(1)/n)
    print(f"  n={n:3d}  |b_n|^(-1/n) = {nstr(r,18)}   r/(pi/2)-1 = {nstr(r/(pi/2)-1,4)}")

# expected correction: |b_n|^{-1/n} = (pi/2) * n^{3/(2n)} * |pat|^{-1/n} (pat has |.|=1)
# extrapolate: r_n*(pi/2)^{-1} = 1 + (3/2) ln n / n + c/n + ...
print("\nfit r_n = (pi/2)(1 + a ln n/n + c/n) on n in [120..180]:")
import numpy as np
ns = list(range(120,181,2))
xs = []; ys=[]
for n in ns:
    r = abs(b[n])**(-mpf(1)/n)/(pi/2)
    xs.append([float(log(mpf(n)) / n), float(1/mpf(n))]); ys.append(float(r-1))
A = np.array(xs); Y=np.array(ys)
coef, res, *_ = np.linalg.lstsq(A, Y, rcond=None)
print("  a = %.6f (expect 1.5),  c = %.6f,  residual rms = %.2e" % (coef[0], coef[1], np.std(Y-A@coef)))
print("elapsed %.1fs"%(time.time()-t0))
