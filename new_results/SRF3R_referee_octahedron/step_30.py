import time
from mpmath import mp, mpf, mpc, log, pi, exp, fabs
t0=time.time()
mp.dps = 50
z0 = (-log(2)+pi*1j)/2
R = abs(z0)
def Q(t,s,u):
    return (1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + u**(2*t)
def dQ(t,s,u):
    return (2*(1-u)**(2*t)*(s**(2*t)*log(s)+(1-s)**(2*t)*log(1-s))
            + 2*log(1-u)*(1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + 2*log(u)*u**(2*t))
def zero_corner(w,U):
    s = exp(-w*U); u = exp(-(1-w)*U)
    # start near the fiber root: use analytic corner approx z0-side root guess
    # solve fiber eq first for a start
    def F(z): return exp(-2*z*w)+exp(-2*z*(1-w))+1
    def dF(z): return -2*w*exp(-2*z*w)-2*(1-w)*exp(-2*z*(1-w))
    z = z0 if w>0.02 else mpc(z0)  # start
    # better: continuation-free Newton from a few starts
    best=None
    for zst in [z0, z0*2, mpc(log(2), pi)]:
        z = zst
        for _ in range(150):
            try: zn = z - F(z)/dF(z)
            except (OverflowError, ZeroDivisionError): break
            if fabs(zn-z) < mpf(10)**(-38): z=zn; break
            z = zn
        if fabs(F(z))<mpf(10)**(-15) and (best is None or fabs(z)<fabs(best)): best=z
    t = best/U
    for _ in range(300):
        tn = t - Q(t,s,u)/dQ(t,s,u)
        if fabs(tn-t) < mpf(10)**(-38): return tn
        t = tn
    return t
print("== actual-Q corner scan: min over w of |tau*| U / R")
for U in [3,4,5,6,8,10,15]:
    rowmin = mpf(9); roww=None
    for i in range(1,40):
        w = mpf(i)/40
        t = zero_corner(w, mpf(U))
        v = fabs(t)*U/R
        if v < rowmin: rowmin, roww = v, w
    print(f"  U={U:3}: min_w |tau*|U/R = {mp.nstr(rowmin,10)} at w={mp.nstr(roww,4)}  [{time.time()-t0:.0f}s]")
print("(fiber-limit values: 2 at w=1/2 down to ~1.0 at edges; actual Q should approach these as U grows)")
