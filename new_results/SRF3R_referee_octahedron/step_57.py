import time
from mpmath import mp, mpf, sqrt, log10, sinh, cosh, tanh, pi, exp, log1p, expm1
mp.dps = 55
t0=time.time()

def ts1d(f, h, N):
    s = mpf(0)
    for j in range(-N, N+1):
        t = mpf(j)*h
        u = (pi/2)*sinh(t)
        x = tanh(u)
        if x == 1 and j > 0: continue
        if x == -1 and j < 0: continue
        w = (pi/2)*cosh(t)/(cosh(u)**2)
        s += w*f((x+1)/2)
    return s*h/2

def onemvPow(v, p):
    # 1 - v^p accurately for v in [0,1]
    if v == 1: return mpf(0)
    return -expm1(p*log1p(v-1)) if v > mpf('0.5') else 1-v**p

def S_graph_ts(p, h, N):
    p = mpf(p)
    def integrand(x, v):
        # z^p = (1-x^p)(1-v^p)  [factored, no cross-cancellation]
        zp = onemvPow(x, p)*onemvPow(v, p)
        z = zp**(1/p)
        zx = -x**(p-1) * z**(1-p); zy = -(ymax_v(x)*v)**(p-1) * z**(1-p)
        return sqrt(1 + zx**2 + zy**2)
    def ymax_v(x):
        return onemvPow(x, p)**(1/p)
    def inner(x):
        ym = ymax_v(x)
        return ym * ts1d(lambda v: integrand(x, v), h, N)
    return 8 * ts1d(inner, h, N)

ST = mpf('6.997417989197152537546960063980097361868')
for h, N in [(mpf('0.2'), 30), (mpf('0.1'), 60), (mpf('0.07'), 90)]:
    v = S_graph_ts(1.01, h, N)
    print(f"h={float(h)}: {mp.nstr(v,48)}")
    print(f"   digits={float(-log10(abs(v-ST)/abs(ST))):.1f}  [{time.time()-t0:.0f}s]")
