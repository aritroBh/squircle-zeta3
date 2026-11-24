import time
from mpmath import mp, mpf, mpc, log, pi, exp, fabs
t0=time.time()
mp.dps = 30
z0 = (-log(2)+pi*1j)/2
R = abs(z0)
def Q(t,s,u):
    return (1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + u**(2*t)
def dQ(t,s,u):
    return (2*(1-u)**(2*t)*(s**(2*t)*log(s)+(1-s)**(2*t)*log(1-s))
            + 2*log(1-u)*(1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + 2*log(u)*u**(2*t))
def fibroot(w):
    # track min-|z| root of fiber eq by continuation in w from 0.5 (root log2+i pi there)
    def F(z): return exp(-2*z*w)+exp(-2*z*(1-w))+1
    def dF(z): return -2*w*exp(-2*z*w)-2*(1-w)*exp(-2*z*(1-w))
    z = mpc(log(2), pi)  # w=0.5 root
    wcur = mpf('0.5'); steps = 40
    dw = (w-wcur)/steps
    for _ in range(steps):
        wcur += dw
        for _ in range(30):
            z = z - (exp(-2*z*wcur)+exp(-2*z*(1-wcur))+1)/(-2*wcur*exp(-2*z*wcur)-2*(1-wcur)*exp(-2*z*(1-wcur)))
    return z
def zero_corner(w,U):
    s = exp(-w*U); u = exp(-(1-w)*U)
    t = fibroot(w)/U
    for _ in range(80):
        tn = t - Q(t,s,u)/dQ(t,s,u)
        if fabs(tn-t) < mpf(10)**(-25): return tn
        t = tn
    return t
print("== actual-Q corner scan: min over w of |tau*| U / R")
for U in [3,5,8,12]:
    rowmin = mpf(9); roww=None
    for i in range(2,20):
        w = mpf(i)/20
        t = zero_corner(w, mpf(U))
        v = fabs(t)*U/R
        if v < rowmin: rowmin, roww = v, w
    print(f"  U={U:3}: min_w |tau*|U/R = {mp.nstr(rowmin,8)} at w={mp.nstr(roww,4)}  [{time.time()-t0:.0f}s]")
