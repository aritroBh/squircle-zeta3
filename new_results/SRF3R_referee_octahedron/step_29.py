import time
from mpmath import mp, mpf, mpc, log, pi, exp, fabs
t0=time.time()
mp.dps = 50
z0 = (-log(2)+pi*1j)/2
R = abs(z0)
def F(z,w): return exp(-2*z*w)+exp(-2*z*(1-w))+1
def dF(z,w): return -2*w*exp(-2*z*w)-2*(1-w)*exp(-2*z*(1-w))
def roots_fiber(w):
    found=[]
    for re0 in [mpf(i)/10 for i in range(-30,31,4)]:
        for im0 in [mpf(i)/10 for i in range(1,50,4)]:
            z = mpc(re0, im0); ok=True
            for _ in range(150):
                try:
                    zn = z - F(z,w)/dF(z,w)
                except (OverflowError, ZeroDivisionError):
                    ok=False; break
                if fabs(zn) > 4*R or fabs(zn-z) > 10*R: ok=False; break
                if fabs(zn-z) < mpf(10)**(-38): z=zn; break
                z = zn
            if ok and fabs(F(z,w)) < mpf(10)**(-18) and fabs(z) < 2.6*R:
                if not any(fabs(z-g)<mpf(10)**(-14) for g in found): found.append(z)
    return found
res=[]
for wi in [0.5,0.4,0.3,0.25,0.2,0.15,0.1,0.07,0.05,0.03,0.02,0.01,0.005]:
    w = mpf(wi); rr = roots_fiber(w)
    mn = min(fabs(z) for z in rr)
    res.append((wi, float(mn/R), len(rr)))
    print(f"  w={wi:6.3f}: min|z|/R = {mp.nstr(mn/R,10)}   (#roots {len(rr)})  [{time.time()-t0:.0f}s]")
print("min over scanned w:", min(r[1] for r in res), " (claim: >=1, =1 only at w->0,1; =2 at w=1/2)")
