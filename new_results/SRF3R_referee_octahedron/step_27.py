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
def zero(s,u,m=-1):
    x = -log(u)
    t = (-log(2)+1j*pi*(2*m+1))/(2*x) if m!=-1 else (-log(2)+1j*pi)/(2*x)
    t = (-log(2)+ (1j*pi*(2*m+1)))/(2*x)
    for _ in range(300):
        tn = t - Q(t,s,u)/dQ(t,s,u)
        if fabs(tn-t) < mpf(10)**(-38): return tn
        t = tn
    return t

print("== TAIL: min over s in [0.05,0.95] of |tau*|x/R, fine grids")
for xe in [5,6,7,8,10,12,15,20]:
    u = exp(-mpf(xe)); rowmin = mpf(9); rows=None
    for i in range(5,96):
        s = mpf(i)/100
        for m in [-1,0]:
            t = zero(s,u,m)
            v = fabs(t)*xe/R
            if v < rowmin: rowmin, rows = v, s
    print(f"  x={xe}: min |tau*|x/R = {mp.nstr(rowmin,10)} at s={mp.nstr(rows,4)}  [{time.time()-t0:.0f}s]")
print("asymptotic prediction: 1 + log2/x (s=1/2 correction):", [round(1+0.6931471805599453/xe,6) for xe in [5,6,7,8,10,12,15,20]])
