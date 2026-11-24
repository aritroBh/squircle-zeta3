import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 50
def spow(t, a):
    if t == 0: return mpf(0)
    return (1 if t>0 else -1)*abs(t)**a
def spowd(t, a):
    if t == 0: return mpf(0)
    return a*abs(t)**(a-1)
def h_oct(c1, c2, p, debug=False):
    q = p/(p-1); sq2 = sqrt(2)
    def G(z3, z4):
        z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
        return abs(z1)**q+abs(z2)**q+abs(z3)**q+abs(z4)**q
    def solve(g3, g4):
        z3, z4 = mpf(g3), mpf(g4)
        for _ in range(80):
            z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
            f1 = -(spow(z1,q-1)+spow(z2,q-1))/sq2 + spow(z3,q-1)
            f2 = -(spow(z1,q-1)-spow(z2,q-1))/sq2 + spow(z4,q-1)
            if abs(f1)+abs(f2) < mpf(10)**(-38): return (z3, z4)
            a1, a2 = spowd(z1,q-1), spowd(z2,q-1)
            J11 = (a1+a2)/2 + spowd(z3,q-1)
            J22 = (a1+a2)/2 + spowd(z4,q-1)
            J12 = (a1-a2)/2; J21 = (a1-a2)/2
            det = J11*J22 - J12*J21
            if det == 0: return None
            dz3 = (-f1*J22 + f2*J12)/det; dz4 = (-J11*f2 + J21*f1)/det
            mx = max(abs(dz3), abs(dz4))
            if mx > mpf('0.3'): dz3, dz4 = dz3*mpf('0.3')/mx, dz4*mpf('0.3')/mx
            z3 += dz3; z4 += dz4
        return None
    r = None
    for g in [(c1/3,c2/3),(c1/2,c2/2),(0,0),(c1,c2),(-c1,-c2),(c1/4,-c2/4),(c2/2,c1/2),(c1,-c2)]:
        r = solve(*g)
        if r: break
    if not r:
        if debug: print("fallback at", nstr(c1,6), nstr(c2,6))
        # nested bisection fallback
        def dGz3(z3, z4):
            z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
            return -(spow(z1,q-1)+spow(z2,q-1))/sq2 + spow(z3,q-1)
        M = 3*max(abs(c1),abs(c2))
        def min_z3(z4):
            lo, hi = -M, M
            for _ in range(110):
                mid = (lo+hi)/2
                if dGz3(mid, z4) < 0: lo = mid
                else: hi = mid
            return G((lo+hi)/2, z4)
        gr = (sqrt(5)-1)/2
        lo, hi = -M, M
        x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
        f1, f2 = min_z3(x1), min_z3(x2)
        for _ in range(50):
            if f1 < f2:
                hi, x2, f2 = x2, x1, f1
                x1 = hi-gr*(hi-lo); f1 = min_z3(x1)
            else:
                lo, x1, f1 = x1, x2, f2
                x2 = lo+gr*(hi-lo); f2 = min_z3(x2)
        return ((f1+f2)/2)**(1/q)
    z3, z4 = r
    return G(z3, z4)**(1/q)
def P_oct(p):
    return 16*mpmath.quad(lambda th: h_oct(cos(th),sin(th),p), [0, pi/16, pi/8])
t0=time.time()
epss = [mpf('0.02')/mpf('1.6')**j for j in range(7)]
Ps = [P_oct(1/e) for e in epss]
json.dump([nstr(e,30) for e in epss], open('ckpt_oct3_epss.json','w'))
json.dump([nstr(P,50) for P in Ps], open('ckpt_oct3_Ps.json','w'))
Poct = 16*(sqrt(2)-1)
for e,P in zip(epss,Ps):
    print(f"eps={nstr(e,10)} Poct-P={nstr(Poct-P,12)}")
print("time", time.time()-t0)
