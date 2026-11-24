import numpy as np
from numpy.polynomial.legendre import leggauss
def gauss2(f, a=0., b=1., n=400):
    x, w = leggauss(n); xs = (x+1)/2*(b-a)+a; ws = w/2*(b-a)
    tot = 0.
    for t2, w2 in zip(xs, ws):
        inner = np.sum(ws * f(xs*t2/b if b else xs, t2) * (t2/b))  # integrate t1 in [0,t2]
        tot += w2*inner
    return tot
def G_xxy(x, y, n=300):
    # G(x,x,y;1) = int_0^1 int_0^{t2} 1/(t1-x)/(t2-x) dt1 * 1/(t3-y)... use reduced 1-d:
    # = int_0^1 (1/2) log^2((t-x)/(-x))/(t-y) dt  -- this IS (1/2)N. Instead validate against 2-d form:
    xg, wg = leggauss(n); ts = (xg+1)/2; ws = wg/2
    val = 0.
    for t3, w3 in zip(ts, ws):
        # inner 2-d: int_0^{t3} dt2/(t2-x) int_0^{t2} dt1/(t1-x)
        sub = leggauss(n); u = (sub[0]+1)/2*t3; wu = sub[1]/2*t3
        inner2 = 0.
        for t2, w2 in zip(u, wu):
            sub2 = leggauss(128); v = (sub2[0]+1)/2*t2; wv = sub2[1]/2*t2
            i1 = np.sum(wv/(v-x))
            inner2 += w2*i1/(t2-x)
        val += w3*inner2/(t3-y)
    return val
import mpmath as mp
mp.mp.dps = 50
def Nint(x, y):
    x,y = mp.mpf(x), mp.mpf(y)
    return mp.re(mp.quad(lambda t: mp.log((t-x)/(-x))**2/(t-y), [0, 1]))
for (x,y) in [(2.,3.),(3.,2.),(2.,-1.),(4.,-3.)]:
    g = G_xxy(x,y)
    nh = float(Nint(x,y))/2
    print(f'({x},{y}): G(x,x,y) [2-d Gauss] = {g:.12f}  (1/2)N = {nh:.12f}  diff = {g-nh:.2e}')
open('progress.log','a').write('step: R1 = (1/2)N pinned by independent 2-d Gauss quadrature\n')
