# SM5: grid check of the N=5 remainder bound for the soft-max window (H2):
#   |L - tau g1 - tau^3 g3| <= C tau^5 (1+s)^6 e^{-s}  on 0<=s<=pi/(4 tau), exact L (with the b term).
import mpmath as mp
mp.mp.dps = 250
worst = (0, None)
for tau in [mp.mpf(x) for x in ("1", "0.5", "0.25", "0.1", "0.05", "0.02", "0.01", "0.005")]:
    smax = mp.pi/(4*tau)
    N = 400
    for i in range(N+1):
        s = smax*i/N
        x = tau*s; u = x
        b = mp.exp(-2*mp.cos(u)/tau)*(1 + mp.exp(-2*mp.sin(u)/tau))
        L = tau*mp.log(1 + mp.exp(-2*mp.sin(x)/tau) + b)
        g1 = mp.log(1 + mp.exp(-2*s)); g3 = s**3/(3*(1 + mp.exp(2*s)))
        r = abs(L - tau*g1 - tau**3*g3)/(tau**5*(1+s)**6*mp.exp(-s))
        if r > worst[0]: worst = (r, (tau, s))
print("largest ratio actual/(tau^5 (1+s)^6 e^{-s}) =", mp.nstr(worst[0], 6), "at (tau, s) =", [mp.nstr(v, 5) for v in worst[1]])
