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

def zero(s,u,branch=1):
    x = -log(u)
    t = (-log(2)+branch*1j*pi)/(2*x)
    for _ in range(300):
        tn = t - Q(t,s,u)/dQ(t,s,u)
        if fabs(tn-t) < mpf(10)**(-40): return tn
        t = tn
    return t

print("== min over s of |tau*| * x / R  (need >= 1 for contour |tau|=(R-eps)/x to be inside zeros)")
glob_min = mpf(9); arg=None
for ue in [-6,-5,-4,-3,-2.5,-2,-1.7,-1.5,-1.3,-1.2,-1.1,-1.0,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.6935,-0.6932]:
    u = exp(mpf(ue)); x = -log(u)
    rowmin = mpf(9); rows=None
    for i in range(1,200):
        s = mpf(i)/200
        t = zero(s,u)
        v = fabs(t)*x/R
        if v < rowmin: rowmin, rows = v, s
    if rowmin < glob_min: glob_min, arg = rowmin, ue
    print(f"  u=e^{ue:7.4f} (x={mp.nstr(x,6)}): min_s |tau*|x/R = {mp.nstr(rowmin,10)} at s={mp.nstr(rows,4)}  [{time.time()-t0:.0f}s]")
print("GLOBAL min:", mp.nstr(glob_min,10), "at log u =", arg)
print("=> zeros stay strictly outside |tau|=R/x for u<1/2; limiting point u=1/2,s=1/2 gives exactly R/x (boundary of the strip).")
