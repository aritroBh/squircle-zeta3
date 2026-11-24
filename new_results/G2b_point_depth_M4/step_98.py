
import mpmath as mp
mp.mp.dps = 100

def J2_quad(b, c, p, t):
    lam = mp.log(abs(p - c)) if p != c else mp.mpf(0)
    f = lambda u: (mp.log(abs(u - c)) - lam) / (u - b) if u != b else mp.mpf(0)
    pts = sorted(set(x for x in [p, t, b, c] if min(p, t) <= x <= max(p, t)))
    if b == p:
        return mp.quad(f, pts)
    # subtract for endpoint/inner poles
    if b == t:
        Ft = mp.quad(lambda u: (mp.log(abs(u-c)) - lam), [p, t]) / 1  # not needed
    if min(p, t) < b < max(p, t):
        fb = mp.log(abs(b - c)) - lam
        g = lambda u: ((mp.log(abs(u - c)) - lam) - fb) / (u - b) if u != b else mp.mpf(0)
        pv = fb * (mp.log(abs(t - b)) - mp.log(abs(b - p)))
        return mp.quad(g, sorted(set([p, b, t]))) + pv
    return mp.quad(f, pts)

def Gfun(u, b, c):
    return mp.log(abs(u - c)) * mp.log(abs((b - u) / (b - c))) + mp.polylog(2, (u - c) / (b - c)).real

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

for (b, c, p, t) in [(2, 3, mp.mpf(1)/2, 1), (3, -2, mp.mpf(-1), mp.mpf(1)/2),
                     (4, 3, mp.mpf(1)/2, mp.mpf(3)/2), (3, 2, mp.mpf(1)/2, 1),
                     (0, 3, mp.mpf(1)/2, 1)]:
    vq = J2_quad(b, c, p, t)
    vc = J2val(b, c, p, t)
    print(f"J2({b},{c};{p}->{t}): quad {mp.nstr(vq,18)} closed {mp.nstr(vc,18)} diff {mp.nstr(abs(vq-vc),3)}")
