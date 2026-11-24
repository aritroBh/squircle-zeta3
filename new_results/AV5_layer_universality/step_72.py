import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, nstr
mp.dps = 80
def h_lame(theta, p, lam):
    if theta == 0 or theta == pi/2: return mpf(1)
    c, s = cos(theta), sin(theta)
    def F(x):
        return c - (1/lam)*x**(p-1)*(1-x**p)**(1/(lam*p)-1)*s
    lo, hi = mpf(0), 1-mpf(10)**(-60)
    for _ in range(260):
        mid = (lo+hi)/2
        if F(mid) > 0: lo = mid
        else: hi = mid
    x = (lo+hi)/2
    return x*c + (1-x**p)**(1/(lam*p))*s
def P_lame_sup(p, lam):
    return 4*mpmath.quad(lambda th: h_lame(th, p, lam), [0, pi/8, pi/4, 3*pi/8, pi/2])
data = json.load(open('ckpt_lame_p1.json'))
t0=time.time()
for lam in ['2','3','0.5']:
    lamv = mpf(lam)
    epss = [mpf('0.02')/mpf('1.6')**j for j in range(10)]
    Ps = [P_lame_sup(1/e, lamv) for e in epss]
    data[f'lam{lam}_pinf2_epss'] = [nstr(e,30) for e in epss]
    data[f'lam{lam}_pinf2_Ps'] = [nstr(P,75) for P in Ps]
    print(f"lam={lam}: 8-P at eps: {[nstr(8-P,10) for P in Ps[:5]]}")
    print("  time", time.time()-t0)
json.dump(data, open('ckpt_lame_p1.json','w'))
