import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 55
def spow(t, a):
    if t == 0: return mpf(0)
    return (1 if t>0 else -1)*abs(t)**a
def h_oct(c1, c2, p):
    q = p/(p-1); sq2 = sqrt(2)
    def dGz3(z3, z4):
        z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
        return -(spow(z1,q-1)+spow(z2,q-1))/sq2 + spow(z3,q-1)
    def G(z3, z4):
        z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
        return abs(z1)**q+abs(z2)**q+abs(z3)**q+abs(z4)**q
    M = 3*max(abs(c1),abs(c2))
    def min_z3(z4):
        lo, hi = -M, M
        for _ in range(130):
            mid = (lo+hi)/2
            if dGz3(mid, z4) < 0: lo = mid
            else: hi = mid
        return G((lo+hi)/2, z4)
    gr = (sqrt(5)-1)/2
    lo, hi = -M, M
    x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
    f1, f2 = min_z3(x1), min_z3(x2)
    for _ in range(70):
        if f1 < f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi-gr*(hi-lo); f1 = min_z3(x1)
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo+gr*(hi-lo); f2 = min_z3(x2)
    return ((f1+f2)/2)**(1/q)
def P_oct(p):
    return 16*mpmath.quad(lambda th: h_oct(cos(th),sin(th),p), [0, pi/16, pi/8])
t0=time.time()
epss = [mpf('0.02')/mpf('1.6')**j for j in range(7)]
Ps = [P_oct(1/e) for e in epss]
json.dump([nstr(e,30) for e in epss], open('ckpt_oct2_epss.json','w'))
json.dump([nstr(P,50) for P in Ps], open('ckpt_oct2_Ps.json','w'))
Poct = 16*(sqrt(2)-1)
for e,P in zip(epss,Ps):
    print(f"eps={nstr(e,10)} Poct-P={nstr(Poct-P,12)}")
print("time", time.time()-t0)
