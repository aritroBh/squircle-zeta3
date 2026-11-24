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
# fixed beta deep inside fan, tau scaling: delta_h should shrink faster than any power
beta = mpf('0.42')   # 2 tan(0.42) = 0.893 < 1
print("delta_h(beta=0.42, inside fan) vs tau:")
for tau in [mpf('0.02'), mpf('0.01'), mpf('0.005'), mpf('0.002')]:
    dh = h_lame_v(beta, 1+tau, mpf(2)) - h_tent_lam2(beta)
    print(f"  tau={tau}: dh={nstr(dh,10)}  dh/tau={nstr(dh/tau,10)}")
# fixed beta outside fan: should be O(tau) bulk
beta2 = mpf('0.75')
print("delta_h(beta=0.75, outside fan) vs tau:")
for tau in [mpf('0.02'), mpf('0.01'), mpf('0.005'), mpf('0.002')]:
    dh = h_lame_v(beta2, 1+tau, mpf(2)) - h_tent_lam2(beta2)
    print(f"  tau={tau}: dh={nstr(dh,10)}  dh/tau={nstr(dh/tau,10)}")
