
import mpmath as mp, time
mp.mp.dps = 55

def hc(s, th):
    c = mp.cos(th)
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)
def quad_sin_t(th, n):
    nn = mp.mpf(n); sq = mp.sqrt(nn)
    def integrand(t):
        x = 1 - t/sq
        fx = mp.sin(mp.pi*x/2)
        yp = nn*mp.exp((nn-1)*mp.log(fx))*(mp.pi/2)*mp.cos(mp.pi*x/2)
        return hc(yp, th)/sq
    return mp.quad(integrand, [0, 3, 15, sq], maxdegree=8)
def ustar(V):
    return mp.findroot(lambda u: u**2 - mp.log(u) - mp.log(V), (mp.log(V))**mp.mpf('0.5'))

beta = mp.pi**2/8; m = 2
t0=time.time()
print("f=sin(pi x/2): higher contact m=2, beta=pi^2/8")
print("leading: (1-c)sqrt(log n/(2 beta n));  refined: (1-c) u*/sqrt(beta n),  V u* e^{-u*^2}=1, V=2 sqrt(beta n)")
rows = []
for name, th in [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3)]:
    c = mp.cos(th)
    for n in [100, 1000, 10000, 100000]:
        q = quad_sin_t(th, n)
        nn = mp.mpf(n)
        leading = (1-c)*mp.sqrt(mp.log(nn)/(m*beta*nn))
        V = 2*mp.sqrt(beta*nn); us = ustar(V)
        refined = (1-c)*us/mp.sqrt(beta*nn)
        rows.append((name, n, q, leading, q/leading, refined, q/refined))
        print(f" th={name:6s} n={n:6d}: 2-L={mp.nstr(q,14)} lead={mp.nstr(leading,10)} r={mp.nstr(q/leading,7)} | u*={mp.nstr(us,7)} ref={mp.nstr(refined,10)} r={mp.nstr(q/refined,7)}")
print('elapsed', time.time()-t0)
import json
json.dump([[r[0], r[1], mp.nstr(r[2],40), mp.nstr(r[3],40), mp.nstr(r[4],40), mp.nstr(r[5],40), mp.nstr(r[6],40)] for r in rows], open('results_sin.json','w'))
