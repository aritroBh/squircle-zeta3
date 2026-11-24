import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, tan, nstr
mp.dps = 70
def h_lame(theta, p, lam):
    if theta == 0 or theta == pi/2: return mpf(1)
    c, s = cos(theta), sin(theta)
    def F(x):
        return c - (1/lam)*x**(p-1)*(1-x**p)**(1/(lam*p)-1)*s
    lo, hi = mpf(0), 1-mpf(10)**(-50)
    for _ in range(220):
        mid = (lo+hi)/2
        if F(mid) > 0: lo = mid
        else: hi = mid
    x = (lo+hi)/2
    return x*c + (1-x**p)**(1/(lam*p))*s
lam = mpf(2)
t0=time.time()
rows = []
for j in range(1, 24):
    th = pi/2 * j/24
    hsq = cos(th)+sin(th)
    vals = []
    for p in [mpf(100), mpf(300)]:
        d = p*(h_lame(th, p, lam) - hsq)
        vals.append(d)
    rows.append((th, vals))
    print(f"theta={nstr(th,8)} F(100)={nstr(vals[0],16)} F(300)={nstr(vals[1],16)}")
print("time", time.time()-t0)
json.dump([[nstr(th,20)]+[nstr(v,40) for v in vv] for th,vv in rows], open('ckpt_lame_profile.json','w'))
