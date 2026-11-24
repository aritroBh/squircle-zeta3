import time
from mpmath import mp, mpf, quad, sqrt, log10, sinh, cosh, tanh, pi, exp
mp.dps = 50
t0=time.time()

def ts1d(f, h, N):
    """manual tanh-sinh over [0,1], nodes j=-N..N step h"""
    s = mpf(0)
    for j in range(-N, N+1):
        t = mpf(j)*h
        u = (pi/2)*sinh(t)
        x = tanh(u)
        w = (pi/2)*cosh(t)/(cosh(u)**2)
        if x == 1 and j > 0: continue
        fval = f((x+1)/2)  # map [-1,1] -> [0,1]
        s += w*fval
    return s*h/2

def S_graph_ts(p, h, N):
    p = mpf(p)
    def integrand(x, y):
        z = (1 - x**p - y**p)**(1/p)
        zx = -x**(p-1) * z**(1-p); zy = -y**(p-1) * z**(1-p)
        return sqrt(1 + zx**2 + zy**2)
    def inner(x):
        ymax = (1 - x**p)**(1/p)
        return ymax * ts1d(lambda v: integrand(x, ymax*v), h, N)
    return 8 * ts1d(inner, h, N)

ST = mpf('6.997417989197152537546960063980097361868')
for h, N in [(mpf('0.2'), 30), (mpf('0.1'), 60), (mpf('0.05'), 120)]:
    v = S_graph_ts(1.01, h, N)
    print(f"h={float(h)}: {mp.nstr(v,42)}  digits={float(-log10(abs(v-ST)/abs(ST))):.1f}  [{time.time()-t0:.0f}s]")
