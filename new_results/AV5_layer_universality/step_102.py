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
def h_hexinf(c1,c2):
    if c1*c2 >= 0: return max(abs(c1),abs(c2))
    return abs(c1)+abs(c2)
# profile near theta=pi/4: scaled coordinate sigma = (theta-pi/4)/eps
p = mpf(400)
eps = 1/p
print("delta_h/eps near pi/4, p=400 (sigma, dh/eps):")
for sig in [-6,-4,-2,-1,0,1,2,4,6]:
    th = pi/4 + eps*sig
    dh = (h_hex(cos(th),sin(th),p) - h_hexinf(cos(th),sin(th)))/eps
    print(f"  {sig:>4} {nstr(dh,16)}")
# and the global shape: dh*eps vs theta
print("global: eps*(dh) shape... dh(θ)/?:")
for j in range(1,8):
    th = pi/4 * j/4
    dh = h_hex(cos(th),sin(th),p) - h_hexinf(cos(th),sin(th))
    print(f"  th={nstr(th,8)} dh={nstr(dh,10)} dh*p={nstr(dh*p,10)}")
