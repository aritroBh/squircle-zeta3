import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 50
def r_oct(phi, p):
    c, s = cos(phi), sin(phi)
    return (c**p + s**p + mpf(2)**(-p/2)*((c+s)**p + (c-s)**p))**(-1/p)
def red(alpha):
    x = (alpha + pi/8) % (pi/4)
    return abs(x - pi/8)
def h_oct(theta, p):
    def hof(phi):
        return r_oct(phi, p)*cos(red(theta-phi))
    lo, hi = mpf(0), pi/4
    gr = (sqrt(5)-1)/2
    x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
    f1, f2 = hof(x1), hof(x2)
    for _ in range(120):
        if f1 < f2:
            lo, x1, f1 = x1, x2, f2
            x2 = lo+gr*(hi-lo); f2 = hof(x2)
        else:
            hi, x2, f2 = x2, x1, f1
            x1 = hi-gr*(hi-lo); f1 = hof(x1)
    return (f1+f2)/2
def P_oct(p):
    return 16*mpmath.quad(lambda th: h_oct(th,p), [0, pi/16, pi/8])
t0=time.time()
epss = [mpf('0.02')/mpf('1.6')**j for j in range(7)]
Ps = [P_oct(1/e) for e in epss]
json.dump([nstr(e,30) for e in epss], open('ckpt_oct4_epss.json','w'))
json.dump([nstr(P,50) for P in Ps], open('ckpt_oct4_Ps.json','w'))
Poct = 16*(sqrt(2)-1)
for e,P in zip(epss,Ps):
    print(f"eps={nstr(e,10)} Poct-P={nstr(Poct-P,12)}")
print("time", time.time()-t0)
