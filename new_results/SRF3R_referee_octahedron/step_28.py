import time
from mpmath import mp, mpf, mpc, log, pi, exp, fabs
t0=time.time()
mp.dps = 50
z0 = (-log(2)+pi*1j)/2
R = abs(z0)

# --- corner fiber equation F(z,w) = e^{-2zw}+e^{-2z(1-w)}+1 = 0; find ALL roots with |z|<2.5R per w
def F(z,w): return exp(-2*z*w)+exp(-2*z*(1-w))+1
def dF(z,w): return -2*w*exp(-2*z*w)-2*(1-w)*exp(-2*z*(1-w))
def roots_fiber(w):
    found=[]
    for re0 in [mpf(i)/10 for i in range(-32,33,4)]:
        for im0 in [mpf(i)/10 for i in range(1,52,4)]:
            z = mpc(re0, im0)
            for _ in range(200):
                zn = z - F(z,w)/dF(z,w)
                if fabs(zn-z) < mpf(10)**(-40): z=zn; break
                z = zn
            if fabs(F(z,w)) < mpf(10)**(-20) and fabs(z) < 2.6*R:
                if not any(fabs(z-g)<mpf(10)**(-15) for g in found): found.append(z)
    return found
print("== fiber equation: min |z|/R over all roots, per w")
res=[]
for wi in [0.5,0.4,0.3,0.25,0.2,0.15,0.1,0.07,0.05,0.03,0.02,0.01]:
    w = mpf(wi); rr = roots_fiber(w)
    mn = min(fabs(z) for z in rr)
    res.append((wi, mn/R, len(rr)))
    print(f"  w={wi:5.2f}: min|z|/R = {mp.nstr(mn/R,10)}   (#roots found {len(rr)})  [{time.time()-t0:.0f}s]")
print("Claim: >=1 everywhere, ->1 only at w->0,1; =2 at w=1/2.  min over w scanned:", min(r[1] for r in res))
