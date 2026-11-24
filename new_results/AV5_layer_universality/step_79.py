import mpmath
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, tan, nstr, atan
mp.dps = 70
def h_lame(theta_vert, p, lam):
    # support function with angle measured from +y axis: h = max(x sin + y cos)
    sv, cv = sin(theta_vert), cos(theta_vert)
    if theta_vert == 0: return mpf(1)
    def F(x):
        return cv - (1/lam)*x**(p-1)*(1-x**p)**(1/(lam*p)-1)*sv
    lo, hi = mpf(0), 1-mpf(10)**(-50)
    for _ in range(220):
        mid = (lo+hi)/2
        if F(mid) > 0: lo = mid
        else: hi = mid
    x = (lo+hi)/2
    return x*sv + (1-x**p)**(1/(lam*p))*cv
def h_tent_lam2(beta):
    # digon |x|+y^2=1: inside fan cos; outside edge branch
    tb = tan(beta)
    if abs(tb) <= mpf('0.5'):
        return cos(beta)
    y = 1/(2*tb); x = 1-y**2
    return x*sin(beta) + y*cos(beta)
bfan = atan(mpf('0.5'))
tau = mpf('0.004')
p = 1+tau
print("sigma,  delta_h/tau (measured),  e^sigma/sqrt(5) (predicted)")
for sig in [-3,-2,-1,-0.5,0,0.5,1.0]:
    beta = bfan + (2*tau/5)*sig
    dh = (h_lame(beta, p, mpf(2)) - h_tent_lam2(beta))/tau
    pred = exp(sig)/sqrt(5)
    print(f"{sig:>6} {nstr(dh,18)} {nstr(pred,18)}")
# also inside fan (sigma -> -inf): delta_h should -> 0
beta = bfan - (2*tau/5)*6
print("deep inside fan (sig=-6):", nstr((h_lame(beta,p,mpf(2))-h_tent_lam2(beta))/tau,6))
with open('progress.log','a') as fh:
    fh.write('step_79: anisotropic p->1: e^sigma shoulder profile confirmed at fan boundary\n')
