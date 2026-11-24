
import mpmath as mp
mp.mp.dps = 100

def ReLi2(z):
    return mp.polylog(2, z).real

def Gfun(u, b, c):
    """antiderivative of log|u-c|/(u-b), b != c."""
    return mp.log(abs(u - c)) * mp.log(abs((b - u) / (b - c))) + ReLi2((u - c) / (b - c))

def J2val(b, c, p, t):
    """regularized J((b,c); p->t) = int_p^t [log|u-c| - lam]/(u-b) du, lam = log|p-c| (0 if p==c)."""
    lam = mp.log(abs(p - c)) if p != c else mp.mpf(0)
    if b == c:
        F = lambda u: mp.mpf(1)/2 * mp.log(abs(u - b))**2 - lam * mp.log(abs(u - b))
        return F(t) - (F(p) if p != b else mp.mpf(0))
    def F(u):
        if u == b:
            return mp.zeta(2) - mp.log(abs(b - c))**2
        return Gfun(u, b, c) - lam * mp.log(abs(u - b))
    return F(t) - F(p)

def J1val(a, p, t):
    return mp.log(abs(t - a)) - (mp.log(abs(p - a)) if p != a else mp.mpf(0))

# validation of J2 against series (convergent from-0 cases)
def J2_series(b, c, x, N=1500):
    """J((b,c); 0->x) via series: innermost c: -sum (x/c)^m/m; prepend b; prepend..."""
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = -(mp.mpf(1)/c)**m / m
    # prepend b
    t = [mp.mpf(0)] * (N + 1)
    aP = mp.mpf(1); acc = mp.mpf(0)
    for n in range(1, N):
        aP *= b; acc += s[n] * aP
        t[n + 1] = -acc / (aP * b * (n + 1))
    return mp.fsum(t[1:])

ok = True
for (b, c, x) in [(2, 3, 1), (3, -2, 1), (-2, 4, 1), (2, 2, 1)]:
    v1 = J2val(b, c, mp.mpf(0), mp.mpf(x))
    v2 = J2_series(b, c, x)
    d = abs(v1 - v2)
    ok &= d < mp.mpf(10) ** (-40)
    print(f"J2({b},{c};0->{x}): closed {mp.nstr(v1,18)} series {mp.nstr(v2,18)} diff {mp.nstr(d,3)}")
# p != 0 cases validated vs difference of series
for (b, c, p, t) in [(2, 3, mp.mpf(1)/2, 1), (3, -2, mp.mpf(-1), mp.mpf(1)/2), (4, 3, mp.mpf(1)/2, mp.mpf(3)/2)]:
    v1 = J2val(b, c, p, t)
    v2 = J2_series(b, c, t) - J2_series(b, c, p) if abs(b) > abs(t) and abs(b) > abs(p) else None
    if v2 is not None:
        print(f"J2({b},{c};{p}->{t}): closed {mp.nstr(v1,18)} series-diff {mp.nstr(v2,18)} diff {mp.nstr(abs(v1-v2),3)}")
    else:
        print(f"J2({b},{c};{p}->{t}): closed {mp.nstr(v1,18)} (no series check)")
