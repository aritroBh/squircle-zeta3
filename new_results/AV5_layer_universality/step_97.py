import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 70
def spow(t, a):  # sign(t)*|t|^a, safe at 0
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
# validate p=2 ellipse and p=50 values
p2 = mpf(2)
print("p=2 check:", nstr(h_hex(cos(mpf('0.5')),sin(mpf('0.5')),p2),30), nstr(sqrt((2-sin(1))/3),30))
p = mpf(50)
for th in [mpf(0), pi/8, pi/4, 3*pi/8, pi/2]:
    print(nstr(th,10), nstr(h_hex(cos(th),sin(th),p),30))
# expected ~ h(0)=1, h(pi/4)=0.7071(1-tiny), h(pi/2)=1
Pfull = mpmath.quad(lambda th: h_hex(cos(th), sin(th), p), [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 2*pi])
print("P(50) =", nstr(Pfull,40))
print("hexagon limit = 2sqrt2+4 =", nstr(2*sqrt(2)+4,30))
