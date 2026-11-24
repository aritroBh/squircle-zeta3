import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr, ellipe
mp.dps = 70
def h_hex(c1, c2, p):
    q = p/(p-1)
    def F(z):
        return (abs(c1-z)**(q-2)*(c1-z)*-1 - abs(c2-z)**(q-2)*(c2-z) + abs(z)**(q-2)*z)
    # F is monotone increasing; bracket [0, min(c1,c2)] expanded
    lo, hi = -2*max(abs(c1),abs(c2)), 2*max(abs(c1),abs(c2))
    for _ in range(160):
        mid = (lo+hi)/2
        if F(mid) < 0: lo = mid
        else: hi = mid
    z = (lo+hi)/2
    return (abs(c1-z)**q + abs(c2-z)**q + abs(z)**q)**(1/q)
# validate p=2
for th in [mpf(0), mpf('0.5'), mpf('1.0')]:
    exact = sqrt((2-sin(2*th))/3)
    print(nstr(th,6), nstr(h_hex(cos(th),sin(th),mpf(2)),30), nstr(exact,30))
t0=time.time()
Pfull = mpmath.quad(lambda th: h_hex(cos(th), sin(th), mpf(2)), [0, pi/2, pi, 3*pi/2, 2*pi])
print("P full =", nstr(Pfull,40), " ellipse:", nstr(4*ellipe(mpf(2)/3),40), " time", time.time()-t0)
# periodicity check
for th in [mpf('0.1'), mpf('0.3')]:
    print("th:", nstr(h_hex(cos(th),sin(th),mpf(2)),25), " th+pi/3:", nstr(h_hex(cos(th+pi/3),sin(th+pi/3),mpf(2)),25))
