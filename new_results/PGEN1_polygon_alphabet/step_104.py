# step_76: find the trilogarithm relation among the 8 Kummer-Spence points via PSLQ
import mpmath as mp
mp.mp.dps = 100
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(abs(1-x))
def pts(x, y):
    x, y = mp.mpf(x), mp.mpf(y)
    return [x, y, y/x, (1-x)/(1-y), (1-x)/(1-x*y), (1-y)/(1-x*y), x*(1-y)/(1-x*y), y*(1-x)/(1-x*y)]
for (xx,yy) in [('0.3','0.6'), ('0.2','0.5')]:
    vals = [L3(z) for z in pts(xx,yy)] + [mp.zeta(3)]
    rel = mp.pslq(vals, tol=mp.mpf('1e-80'), maxcoeff=1000)
    print(f"({xx},{yy}):", rel)
