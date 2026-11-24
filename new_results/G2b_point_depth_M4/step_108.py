
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 80
from g2b_dbg import s05
C0 = Fr(1,2)

# reproduce the machinery inline
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

def J3_quad(a, b, c, p, q):
    """ground truth via direct quad (mpmath, careful splitting) at dps 60"""
    if p == q: return mp.mpf(0)
    lo, hi = min(p, q), max(p, q)
    inside = sorted(set(x for x in (a, b, c) if lo < x < hi))
    pts = [p] + inside + [q]
    if a == p:
        f = lambda t: J2val(b, c, p, t) / (t - a)
        return mp.quad(f, pts)
    if a == q:
        fq = J2val(b, c, p, q)
        g = lambda t: (J2val(b, c, p, t) - fq) / (t - a) if t != a else mp.mpf(0)
        return mp.quad(g, pts) + fq * mp.log(abs(q - p))
    if lo < a < hi:
        fa = J2val(b, c, p, a)
        g = lambda t: (J2val(b, c, p, t) - fa) / (t - a) if t != a else mp.mpf(0)
        pv = fa * (mp.log(abs(q - a)) - mp.log(abs(a - p)))
        return mp.quad(g, pts) + pv
    f = lambda t: J2val(b, c, p, t) / (t - a)
    return mp.quad(f, pts)

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

def J3_ts(a, b, c, p, q, level=5):
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

# Li5 word terms
w = (0, 0, 0, 0, 2)
e, c = Fr(0), Fr(1, 2)
pts = (e,) + tuple(w) + (c,)
print("points:", pts)
for k in (2, 3):
    for i in range(6 - k):
        gap = w[i:i + k]; ret = w[:i] + w[i + k:]
        ei, ci = pts[i], pts[i + k + 1]
        if k == 2:
            lam = s05.lam2(gap, ci, ei)
            if lam:
                print(f"k=2 i={i}: gap={gap} {ei}->{ci} lam={lam}; ret={ret} (J via series)")
        else:
            lam = s05.lam2(ret, c, e)
            if lam:
                a3, b3, c3 = gap
                v_ts = J3_ts(a3, b3, c3, mp.mpf(ei.numerator)/ei.denominator, mp.mpf(ci.numerator)/ci.denominator)
                v_q = J3_quad(a3, b3, c3, mp.mpf(ei.numerator)/ei.denominator, mp.mpf(ci.numerator)/ci.denominator)
                print(f"k=3 i={i}: gap={gap} {ei}->{ci} lam={lam}: J3 ts={mp.nstr(v_ts,12)} quad={mp.nstr(v_q,12)} diff={mp.nstr(abs(v_ts-v_q),3)}")
