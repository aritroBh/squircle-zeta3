
import mpmath as mp
mp.mp.dps = 80
exec(open("<project>/research/08_theory/G2b_point_depth_M4/step_112.py").read().split("import itertools, time")[0]) if False else None
# rebuild quickly (rely on previous cell defs)
def reli2(z):
    z = mp.mpf(z)
    if z == 0: return mp.mpf(0)
    if z == 1: return mp.pi**2 / 6
    if abs(z) <= mp.mpf('0.55'):
        tot = mp.mpf(0); zk = z; k = 1
        while True:
            t = zk / (k * k); tot += t
            if abs(t) < mp.mpf(10) ** (-mp.mp.dps + 10): break
            zk *= z; k += 1
        return tot
    if z > 1: return -reli2(1 / z) + mp.pi**2 / 3 - mp.log(z)**2 / 2
    if z > 0: return mp.pi**2 / 6 - mp.log(z) * mp.log(1 - z) - reli2(1 - z)
    w = z / (z - 1)
    return -mp.log(1 - z)**2 / 2 - reli2(w)

def Gfun(u, b, c):
    if u == c: return mp.mpf(0)
    if u == b: return mp.zeta(2) - mp.log(abs(b - c))**2
    return mp.log(abs(u - c)) * mp.log(abs((b - u) / (b - c))) + reli2((u - c) / (b - c))

def J2val(b, c, p, t):
    lam = mp.log(abs(p - c)) if p != c else mp.mpf(0)
    if b == c:
        F = lambda u: mp.mpf(1)/2 * mp.log(abs(u - b))**2 - lam * mp.log(abs(u - b))
        return F(t) - (F(p) if p != b else mp.mpf(0))
    def F(u):
        if u == b: return mp.zeta(2) - mp.log(abs(b - c))**2
        return Gfun(u, b, c) - lam * mp.log(abs(u - b))
    return F(t) - F(p)

def ts(f, a, b, level=5):
    if a == b: return mp.mpf(0)
    sgn = 1
    if b < a: a, b = b, a; sgn = -1
    h = mp.mpf(2) ** (-level)
    K = int(mp.mpf('4.6') * 2**level)
    S = mp.mpf(0); pih = mp.pi / 2
    for k in range(-K, K + 1):
        t = k * h
        sh = mp.sinh(t)
        u = mp.tanh(pih * sh)
        x = a + (b - a) * (1 + u) / 2
        dx = (b - a) / 2 * pih * mp.cosh(t) / (mp.cosh(pih * sh) ** 2)
        S += f(x) * dx
    return sgn * h * S

def J3val(a, b, c, p, q, level=5):
    if p == q: return mp.mpf(0)
    if a == p and b == a and c == a:
        return mp.log(abs(q - p))**3 / 6 * (1 if q > p else -1)
    lo, hi = min(p, q), max(p, q)
    inside = sorted(set(x for x in (b, c) if lo < x < hi and x != a))
    pts = [p] + inside + [q]
    tot = mp.mpf(0)
    if a == p:
        for u, v in zip(pts, pts[1:]):
            tot += ts(lambda t: J2val(b, c, p, t) / (t - a), u, v, level)
        return tot
    if a == q:
        fq = J2val(b, c, p, q)
        for u, v in zip(pts, pts[1:]):
            tot += ts(lambda t: (J2val(b, c, p, t) - fq) / (t - a) if t != a else mp.mpf(0), u, v, level)
        return tot + fq * mp.log(abs(q - p))
    if lo < a < hi:
        fa = J2val(b, c, p, a)
        for u, v in zip(pts, pts[1:]):
            tot += ts(lambda t: (J2val(b, c, p, t) - fa) / (t - a) if t != a else mp.mpf(0), u, v, level)
        return tot + fa * (mp.log(abs(q - a)) - mp.log(abs(a - p)))
    for u, v in zip(pts, pts[1:]):
        tot += ts(lambda t: J2val(b, c, p, t) / (t - a), u, v, level)
    return tot

h = mp.mpf(1)/2
j100 = J3val(1, 0, 0, mp.mpf(0), h)
j010 = J3val(0, 1, 0, mp.mpf(0), h)
j001 = J3val(0, 0, 1, mp.mpf(0), h)
j1 = mp.log(1 - h)
j00 = mp.log(h)**2 / 2
print("J(1,0,0) =", mp.nstr(j100, 14))
print("J(0,1,0) =", mp.nstr(j010, 14))
print("J(0,0,1) =", mp.nstr(j001, 14))
print("shuffle sum:", mp.nstr(j100 + j010 + j001, 12), "vs J(1)*J(0,0):", mp.nstr(j1 * j00, 12))
print("diff:", mp.nstr(abs(j100 + j010 + j001 - j1*j00), 3))
