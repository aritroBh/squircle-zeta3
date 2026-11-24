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
        for _ in range(140):
            mid = (lo+hi)/2
            if dGz3(mid, z4) < 0: lo = mid
            else: hi = mid
        return G((lo+hi)/2, z4)
    gr = (sqrt(5)-1)/2
    lo, hi = -M, M
    x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
    f1, f2 = min_z3(x1), min_z3(x2)
    for _ in range(80):
        if f1 < f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi-gr*(hi-lo); f1 = min_z3(x1)
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo+gr*(hi-lo); f2 = min_z3(x2)
    return ((f1+f2)/2)**(1/q)
p = mpf(30)
print("period pi/2 check:", nstr(h_oct(cos(mpf('0.1')),sin(mpf('0.1')),p),18), nstr(h_oct(cos(mpf('0.1')+pi/2),sin(mpf('0.1')+pi/2),p),18))
print("period pi/4 check:", nstr(h_oct(cos(mpf('0.1')),sin(mpf('0.1')),p),18), nstr(h_oct(cos(mpf('0.1')+pi/4),sin(mpf('0.1')+pi/4),p),18))
def P_oct(p):
    return 8*mpmath.quad(lambda th: h_oct(cos(th),sin(th),p), [0, pi/8, pi/4])
t0=time.time()
print("P(2):", nstr(P_oct(mpf(2)),25), nstr(pi*sqrt(2),25), "time", time.time()-t0)
