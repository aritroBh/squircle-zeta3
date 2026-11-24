
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 80
from g2b_dbg import s05
C0 = Fr(1,2)

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

def J3_ts(a, b, c, p, q, level=5):
    if p == q: return mp.mpf(0)
    # analytic case: all three letters == base point (tangential regularization)
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

# retest the Li5 word
w = (0, 0, 0, 0, 2)
e, c = Fr(0), Fr(1, 2)
pts = (e,) + tuple(w) + (c,)
tot = mp.mpf(0)
for k in (2, 3):
    for i in range(6 - k):
        gap = w[i:i + k]; ret = w[:i] + w[i + k:]
        ei, ci = pts[i], pts[i + k + 1]
        if k == 2:
            lam = s05.lam2(gap, ci, ei)
            if lam:
                # from-0 retained via HPL series
                NSER = 800
                aa = ret[-1]
                s = [mp.mpf(0)] * (NSER + 1)
                for m in range(1, NSER + 1):
                    s[m] = mp.mpf(1)/m if aa == 1 else mp.mpf((-1)**(m+1))/m
                for a4 in reversed(ret[:-1]):
                    t = [mp.mpf(0)] * (NSER + 1)
                    if a4 == 0:
                        for m in range(1, NSER + 1): t[m] = s[m]/m
                    else:
                        acc = mp.mpf(0)
                        for m in range(0, NSER):
                            acc = acc * a4 + s[m]; t[m+1] = acc/(m+1)
                    s = t
                Hv = mp.fsum(s[m] * (mp.mpf(1)/2)**m for m in range(1, NSER + 1))
                Jv = ((-1) ** sum(1 for a4 in ret if a4 == 1)) * Hv
                tot += lam * Jv
                print(f"k=2 i={i}: lam={lam} J(ret={ret};0->1/2)={mp.nstr(Jv,10)}")
        else:
            lam = s05.lam2(ret, c, e)
            if lam:
                v = J3_ts(gap[0], gap[1], gap[2], mp.mpf(ei.numerator)/ei.denominator, mp.mpf(ci.numerator)/ci.denominator)
                tot -= lam * v
                print(f"k=3 i={i}: gap={gap} {ei}->{ci} lam={lam} J3={mp.nstr(v,10)}")
print("delta_num(Li5(1/4) word) =", mp.nstr(tot, 10))
