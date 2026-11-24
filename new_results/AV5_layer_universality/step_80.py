import mpmath
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, tan, nstr, atan
mp.dps = 70
def h_lame_v(beta, p, lam):
    sv, cv = sin(beta), cos(beta)
    if beta == 0: return mpf(1)
    def F(x):
        return sv - (1/lam)*x**(p-1)*(1-x**p)**(1/(lam*p)-1)*cv
    lo, hi = mpf(0), 1-mpf(10)**(-50)
    for _ in range(220):
        mid = (lo+hi)/2
        if F(mid) > 0: lo = mid
        else: hi = mid
    x = (lo+hi)/2
    return x*sv + (1-x**p)**(1/(lam*p))*cv
def h_tent_lam2(beta):
    tb = tan(beta)
    if abs(tb) <= mpf('0.5'):
        return cos(beta)
    y = 1/(2*tb); x = 1-y**2
    return x*sin(beta) + y*cos(beta)
# sanity: at p=1, h_lame_v should equal h_tent_lam2
for beta in [mpf('0.2'), mpf('0.4636476090008061'), mpf('0.7')]:
    print(nstr(beta,8), nstr(h_lame_v(beta, mpf(1), mpf(2)),25), nstr(h_tent_lam2(beta),25))
bfan = atan(mpf('0.5'))
tau = mpf('0.004')
p = 1+tau
print("\nsigma,  delta_h/tau (measured),  e^sigma/sqrt(5) (predicted)")
for sig in [-4,-3,-2,-1,-0.5,0,0.5]:
    beta = bfan + (2*tau/5)*sig
    dh = (h_lame_v(beta, p, mpf(2)) - h_tent_lam2(beta))/tau
    pred = exp(sig)/sqrt(5)
    print(f"{sig:>6} {nstr(dh,18)} {nstr(pred,18)}")
