
import mpmath as mp, time
mp.mp.dps = 130

def reli2(z):
    z = mp.mpf(z)
    if z == 0: return mp.mpf(0)
    if abs(z) <= mp.mpf('0.55'):
        tot = mp.mpf(0); zk = z; k = 1
        while True:
            t = zk / (k * k)
            tot += t
            if abs(t) < mp.mpf(10) ** (-mp.mp.dps + 10): break
            zk *= z; k += 1
        return tot
    if z > 1:
        return -reli2(1 / z) + mp.pi**2 / 3 - mp.log(z)**2 / 2
    if z > 0:
        return mp.pi**2 / 6 - mp.log(z) * mp.log(1 - z) - reli2(1 - z)
    w = z / (z - 1)
    return -mp.log(1 - z)**2 / 2 - reli2(w)

def ts(f, a, b, level):
    if a == b: return mp.mpf(0)
    sgn = 1
    if b < a: a, b = b, a; sgn = -1
    h = mp.mpf(2) ** (-level)
    K = int(mp.mpf('4.6') * 2**level)
    S = mp.mpf(0)
    pih = mp.pi / 2
    for k in range(-K, K + 1):
        t = k * h
        sh = mp.sinh(t)
        ch = mp.cosh(mp.pi/2 * sh) if abs(sh) < 200 else None
        u = mp.tanh(pih * sh)
        x = a + (b - a) * (1 + u) / 2
        dx = (b - a) / 2 * pih * mp.cosh(t) / (mp.cosh(pih * sh) ** 2)
        S += f(x) * dx
    return sgn * h * S

def Gfun(u, b, c):
    return mp.log(abs(u - c)) * mp.log(abs((b - u) / (b - c))) + reli2((u - c) / (b - c))

def J2val(b, c, p, t):
    lam = mp.log(abs(p - c)) if p != c else mp.mpf(0)
    if b == c:
        F = lambda u: mp.mpf(1)/2 * mp.log(abs(u - b))**2 - lam * mp.log(abs(u - b))
        return F(t) - (F(p) if p != b else mp.mpf(0))
    def F(u):
        if u == b:
            return mp.zeta(2) - mp.log(abs(b - c))**2
        return Gfun(u, b, c) - lam * mp.log(abs(u - b))
    return F(t) - F(p)

# accuracy test: J3(2,3,4; 0->1) = ∫₀¹ J2(3,4;0,t)/(t-2) dt, vs mpmath quad at high precision
def J3val(a, b, c, p, q, level=8):
    if p == q: return mp.mpf(0)
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

mp.mp.dps = 200
ref = mp.quad(lambda t: J2val(3, 4, mp.mpf(0), t) / (t - 2), [0, 1])
mp.mp.dps = 130
for L in (7, 8, 9):
    t0 = time.time()
    v = J3val(2, 3, 4, mp.mpf(0), mp.mpf(1), L)
    print(f"level {L}: diff {mp.nstr(abs(v-ref),3)}  ({time.time()-t0:.1f}s)")
