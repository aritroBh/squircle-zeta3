import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 70
def spow(t, a):
    if t == 0: return mpf(0)
    return (1 if t>0 else -1)*abs(t)**a
def h_hex(c1, c2, p):
    q = p/(p-1)
    def F(z):
        return -spow(c1-z,q-1) - spow(c2-z,q-1) + spow(z,q-1)
    lo, hi = -2*max(abs(c1),abs(c2)), 2*max(abs(c1),abs(c2))
    for _ in range(160):
        mid = (lo+hi)/2
        if F(mid) < 0: lo = mid
        else: hi = mid
    z = (lo+hi)/2
    return (abs(c1-z)**q + abs(c2-z)**q + abs(z)**q)**(1/q)
def P_hex(p):
    return 4*mpmath.quad(lambda th: h_hex(cos(th), sin(th), p), [0, pi/8, pi/4, 3*pi/8, pi/2])
t0=time.time()
epss = [mpf('0.02')/mpf('1.55')**j for j in range(8)]
Ps = [P_hex(1/e) for e in epss]
json.dump([nstr(e,30) for e in epss], open('ckpt_hex2_epss.json','w'))
json.dump([nstr(P,65) for P in Ps], open('ckpt_hex2_Ps.json','w'))
Phex = 2*sqrt(2)+4
for e,P in zip(epss,Ps):
    print(f"eps={nstr(e,10)} Phex-P={nstr(Phex-P,15)}")
print("time", time.time()-t0)
