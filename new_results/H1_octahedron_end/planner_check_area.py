# Planner check of H1: direct surface area of |x|^p+|y|^p+|z|^p=1 vs S(1+d) = 4 sqrt3 [1 + d + (1-pi^2/9)d^2 + (1+pi^2/9-16 z3/9)d^3 + (1-2pi^2/9+32 z3/9-16 pi^4/405)d^4 + ...]
import mpmath as mp
mp.mp.dps = 22
def S(p):
    # graph over the first-octant face region D = {x^p + y^p <= 1}, with x = u, y = (1-u^p)^(1/p) v
    def f(u, v):
        Y = (1 - u**p)**(1/p); y = Y*v
        w = 1 - u**p - y**p
        if w <= 0: return mp.mpf(0)
        zx = u**(p-1)*w**(1/p-1); zy = y**(p-1)*w**(1/p-1)
        return mp.sqrt(1 + zx*zx + zy*zy)*Y
    return 8*mp.quad(f, [0, 1], [0, 1], maxdegree=8)
z3, pi = mp.zeta(3), mp.pi
def ser(d, K):
    c = [1, 1, 1 - pi**2/9, 1 + pi**2/9 - 16*z3/9, 1 - 2*pi**2/9 + 32*z3/9 - 16*pi**4/405]
    return 4*mp.sqrt(3)*sum(c[k]*d**k for k in range(K+1))
for d in (mp.mpf('0.04'), mp.mpf('0.02'), mp.mpf('0.01')):
    s = S(1 + d)
    print("d", d, "direct", mp.nstr(s, 16), " err vs series to d^3:", mp.nstr(s - ser(d, 3), 4), " to d^4:", mp.nstr(s - ser(d, 4), 4))
