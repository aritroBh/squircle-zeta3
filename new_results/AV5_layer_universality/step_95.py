import mpmath
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 70
def h_hex(c1, c2, p):
    q = p/(p-1)
    def F(z):
        return (-abs(c1-z)**(q-2)*(c1-z) - abs(c2-z)**(q-2)*(c2-z) + abs(z)**(q-2)*z)
    lo, hi = -2*max(abs(c1),abs(c2)), 2*max(abs(c1),abs(c2))
    for _ in range(160):
        mid = (lo+hi)/2
        if F(mid) < 0: lo = mid
        else: hi = mid
    z = (lo+hi)/2
    return (abs(c1-z)**q + abs(c2-z)**q + abs(z)**q)**(1/q)
p = mpf(50)
for th in [mpf(0), pi/8, pi/4, 3*pi/8, pi/2]:
    print(nstr(th,10), nstr(h_hex(cos(th),sin(th),p),30))
# expected: h(0)~1, h(pi/4)~0.7071, h(pi/2)~1
# hexagon support function for comparison
def h_hexinf(c1,c2):
    return max([c1*1+c2*0, c1*0+c2*1, c1*(-1)+c2*1, c1*(-1)+c2*0, c1*0+c2*(-1), c1*1+c2*(-1)])
for th in [mpf(0), pi/8, pi/4, 3*pi/8, pi/2]:
    print("inf:", nstr(th,10), nstr(h_hexinf(cos(th),sin(th)),30))
