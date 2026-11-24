import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, log, log10, nstr, cos
import numpy as np
mp.dps = 150
f450 = [mpf(v) for v in json.load(open('step2_num450_N180.json'))['f450']]
sq2 = sqrt(2)
f1 = [f450[n] - sq2*((-1)**n) for n in range(len(f450))]
b  = [f450[n]/factorial(n) for n in range(len(f450))]

# ---- extrapolated root test ----
ns = list(range(120,181,2))
xs = []; ys=[]
for n in ns:
    r = abs(b[n])**(-mpf(1)/n)/(pi/2)
    xs.append([float(log(mpf(n)) / n), float(1/mpf(n))]); ys.append(float(r-1))
A = np.array(xs); Y=np.array(ys)
coef, res, *_ = np.linalg.lstsq(A, Y, rcond=None)
print("root-test fit r_n/(pi/2) = 1 + a ln n/n + c/n : a=%.5f (theory 1.5), c=%.5f, rms=%.2e"%(coef[0],coef[1],np.std(Y-A@coef)))

# ---- Domb-Sykes with period 4:  r_n = (b_n/b_{n+4})^{1/4} = R (1 + a/n + c/n^2 + d/n^3) ----
# b_n ~ pat(n) (2/pi)^n n^{-3/2}; pat(n)=pat(n+4), so ratio isolates radius R=pi/2 and exponent 3/2
def DS(n): return (b[n]/b[n+4])**mpf('0.25')
print("\nDomb-Sykes r_n=(b_n/b_{n+4})^{1/4}, values:")
for n in [100,120,140,160,170]:
    print(f"  n={n}: {nstr(DS(n),20)}   /(pi/2)-1 = {nstr(DS(n)/(pi/2)-1,4)}")

def dsfit(ns, p):
    # r_n = R*(1 + a1/n + a2/n^2 + ... + ap/n^p)
    M = []
    y = []
    for n in ns:
        M.append([mpf(n)**-k for k in range(1,p+1)])
        y.append(DS(n))
    Mm = mp.matrix(M); ym = mp.matrix(y)
    # least squares via normal equations at mp precision
    G = Mm.T*Mm; h = Mm.T*ym
    x = (G**-1)*h
    # R from constraint r_n/(1+...) : better fit R explicitly: r_n = R + R*a1/n... linear in (R,Ra_i)
    M2 = []
    for n in ns:
        M2.append([mpf(1)]+[mpf(n)**-k for k in range(1,p+1)])
    Mm = mp.matrix(M2); 
    G = Mm.T*Mm; h = Mm.T*ym
    x = (G**-1)*h
    return x  # x[0]=R, x[1]=R*a1,...
for ns,p in [([150,155,160,165,170,175],3),([140,150,160,170,175],3),([160,165,170,175],2),([120,130,140,150,160,170],4)]:
    x = dsfit(ns,p)
    R = x[0]; a1 = x[1]/R if p>=1 else None
    print(f"  DS fit ns={ns} p={p}: R={nstr(R,22)}  err_vs_pi/2={nstr(R-pi/2,3)}  a1={nstr(a1,8)} (theory 1.5)")
