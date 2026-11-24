
import sys, time
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 80

def Jval(w, p, q, cache):
    w = tuple(w)
    if not w: return mp.mpf(1)
    if len(w) == 1:
        a = w[0]
        return mp.log(abs(q - a)) - (mp.log(abs(p - a)) if p != a else mp.mpf(0))
    if p == q: return mp.mpf(0)
    a = w[0]
    lo, hi = min(p, q), max(p, q)
    inside = sorted(set(b for b in w[1:] if lo < b < hi))
    if a == p:
        f = lambda t: Jval(w[1:], p, t, cache) / (t - a)
        return mp.quad(f, [p] + inside + [q])
    if a == q:
        fq = Jval(w[1:], p, q, cache)
        g = lambda t: (Jval(w[1:], p, t, cache) - fq) / (t - a)
        return mp.quad(g, [p] + inside + [q]) + fq * mp.log(abs(q - p))
    if lo < a < hi:
        fa = Jval(w[1:], p, a, cache)
        def g(t):
            return (Jval(w[1:], p, t, cache) - fa) / (t - a)
        pv = fa * (mp.log(abs(q - a)) - mp.log(abs(a - p)))
        return mp.quad(g, [p, a] + [b for b in inside if b != a] + [q]) + pv
    f = lambda t: Jval(w[1:], p, t, cache) / (t - a)
    return mp.quad(f, [p] + inside + [q])

# validation 1: J((0,0,1); 0->1/2) vs hpl series (= -H_{0,0,1}(1/2)? sign: J = (-1)^1 H? H = (-1)^{#1} J)
t0 = time.time()
v1 = Jval((0, 0, 1), mp.mpf(0), mp.mpf(1)/2, {})
# hpl series
def hpl(w, x, N):
    a = w[-1]
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
    for a in reversed(w[:-1]):
        t = [mp.mpf(0)] * (N + 1)
        if a == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N):
                acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[m] * x ** m for m in range(1, N + 1))
h1 = -hpl((0, 0, 1), mp.mpf(1)/2, 500)   # J = (-1)^{#1=1} H -> J = -H
print("J((0,0,1);0->1/2) quad:", mp.nstr(v1, 25))
print("             series:", mp.nstr(h1, 25), "  err:", mp.nstr(abs(v1-h1), 3), f"({time.time()-t0:.2f}s)")

# validation 2: Li32 word: J((0,0,2,0,3); 0->1) should be -Li32(1/2, 1/3)? 1/x=2 -> x=1/2, 1/(xy)=3 -> xy=1/3 -> y=2/3
t0 = time.time()
v2 = Jval((0, 0, 2, 0, 3), mp.mpf(0), mp.mpf(1), {})
def li32(x, y, digits=80):
    x = mp.mpf(x); y = mp.mpf(y)
    Nl = int((digits + 20) * mp.log(10) / -mp.log(abs(x*y))) + 50
    S = mp.mpf(0); ym = y; xn = x; tot = mp.mpf(0)
    for n in range(2, Nl + 1):
        S += ym / ((n-1)*(n-1)); ym *= y; xn *= x
        tot += xn * S / (n*n*n)
    return tot
l2 = li32(mp.mpf(1)/2, mp.mpf(2)/3)
print("J((0,0,2,0,3);0->1) quad:", mp.nstr(v2, 25))
print("        -Li32(1/2,2/3):", mp.nstr(-l2, 25), "  err:", mp.nstr(abs(v2+l2), 3), f"({time.time()-t0:.2f}s)")

# validation 3: a PV case: J((0,-1,0); -1->1) -- check via shuffle relation? just time it
t0 = time.time()
v3 = Jval((0, -1, 0), mp.mpf(-1), mp.mpf(1), {})
print("J((0,-1,0);-1->1) PV:", mp.nstr(v3, 25), f"({time.time()-t0:.2f}s)")
