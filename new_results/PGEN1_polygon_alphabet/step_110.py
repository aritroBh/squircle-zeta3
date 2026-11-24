# step_82: hunt the exact Kummer-Spence relation; test on L3
import mpmath as mp
mp.mp.dps = 100
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(abs(1-x))
def P(x, y):
    """the 8 Kummer-Spence points"""
    return [x, y, y/x, (1-x)/(1-y), (1-x)/(1-x*y), (1-y)/(1-x*y), x*(1-y)/(1-x*y), y*(1-x)/(1-x*y)]
# candidate sign patterns from literature (Goncharov math/0103059 Prop 2.3):
# +x -y +y/x +(1-x)/(1-y) -(1-x)/(1-xy) +(1-y)/(1-xy) -(x(1-y)/(1-xy)) +(y(1-x)/(1-xy))
pat = [1,-1,1,1,-1,1,-1,1]
for xx, yy in [('0.3','0.6'), ('0.2','0.7'), ('0.45','0.85')]:
    x, y = mp.mpf(xx), mp.mpf(yy)
    pts = P(x,y)
    val = mp.fsum([c*L3(z) for c,z in zip(pat, pts)])
    print(f"KS({xx},{yy}) with Goncharov pattern =", mp.nstr(val, 8))
# solve for the true sign pattern: pslq on the 8 values for a non-degenerate pair
x, y = mp.mpf('0.31'), mp.mpf('0.72')
vals = [L3(z) for z in P(x,y)]
rel = mp.pslq(vals, tol=mp.mpf('1e-80'), maxcoeff=100)
print("pslq on 8 KS points (0.31,0.72):", rel)
x, y = mp.mpf('0.17'), mp.mpf('0.53')
vals = [L3(z) for z in P(x,y)]
rel = mp.pslq(vals, tol=mp.mpf('1e-80'), maxcoeff=100)
print("pslq on 8 KS points (0.17,0.53):", rel)
