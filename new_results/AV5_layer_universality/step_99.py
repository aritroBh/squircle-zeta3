import mpmath
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
p = mpf(50)
for th in [mpf('0.2'), mpf('0.6'), mpf('1.0'), mpf('1.3')]:
    h0 = h_hex(cos(th),sin(th),p)
    hA = h_hex(cos(th+pi/2),sin(th+pi/2),p)
    hB = h_hex(cos(pi/2-th),sin(pi/2-th),p)
    hC = h_hex(cos(pi-th),sin(pi-th),p)
    print(f"th={nstr(th,8)}: {nstr(h0,16)}  +90deg:{nstr(hA,16)}  swap:{nstr(hB,16)}  pi-th:{nstr(hC,16)}")
