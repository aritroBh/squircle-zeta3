
import mpmath as mp
mp.mp.dps = 100

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

def J_series(w, x, N=1500):
    """J(word; 0->x) series; letters nonzero or 0, rightmost nonzero, |letters| >= |x|."""
    w = list(w)
    a = w[-1]
    s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1):
        s[m] = -(mp.mpf(1)/a)**m / m
    for letter in reversed(w[:-1]):
        if letter == 0:
            s = [mp.mpf(0)] + [s[m]/m for m in range(1, N)] + [mp.mpf(0)]
        else:
            t = [mp.mpf(0)] * (N + 1)
            aP = mp.mpf(1); acc = mp.mpf(0)
            for n in range(1, N):
                aP *= letter; acc += s[n] * aP
                t[n + 1] = -acc / (aP * letter * (n + 1))
            s = t
    return mp.fsum(s[1:])

def J_antipode(w, p, t, N=1200):
    """J((a1..an); p->t) = sum_i (-1)^{n-i} J(a1..ai; 0->t) J(p-a_n,...,p-a_{i+1}; 0->p)."""
    n = len(w)
    tot = mp.mpf(0)
    for i in range(n + 1):
        left = J_series(w[:i], t, N) if i else mp.mpf(1)
        rw = [p - a for a in reversed(w[i:])]
        right = J_series(rw, p, N) if rw else mp.mpf(1)
        tot += (-1) ** (n - i) * left * right
    return tot

for (b, c, p, t) in [(2, 3, mp.mpf(1)/2, 1), (3, -2, mp.mpf(-1), mp.mpf(1)/2), (4, 3, mp.mpf(1)/2, mp.mpf(3)/2)]:
    v1 = J2val(b, c, p, t)
    v2 = J_antipode((b, c), p, t)
    print(f"J2({b},{c};{p}->{t}): closed {mp.nstr(v1,18)} antipode {mp.nstr(v2,18)} diff {mp.nstr(abs(v1-v2),3)}")
