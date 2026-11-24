import mpmath as mp

def setdps(n):
    mp.mp.dps = n

def Li_ab(a, b, x, y, terms=400):
    """Li_{a,b}(x,y) = sum_{0<k1<k2} x^k1 y^k2/(k1^a k2^b); series."""
    x = mp.mpf(x); y = mp.mpf(y)
    s = mp.mpf(0)
    for k2 in range(2, terms):
        inner = mp.mpf(0)
        xp = x
        for k1 in range(1, k2):
            inner += xp / mp.power(k1, a)
            xp *= x
        s += y**k2 / mp.power(k2, b) * inner
    return s

def Li_n(n, z, terms=2000):
    z = mp.mpf(z)
    s = mp.mpf(0); zp = z
    for k in range(1, terms):
        s += zp/mp.power(k, n)
        zp *= z
    return s

def G_int(args, z):
    """G(a1,...,an; z) = int_0^z dt1/(t1-a1) int_0^{t1} ... ; outermost = first arg."""
    args = [mp.mpf(a) for a in args]
    z = mp.mpf(z)
    def rec(idx, lo):
        if idx == len(args):
            return mp.mpf(1)
        a = args[idx]
        return mp.quad(lambda t: rec(idx+1, t)/(t-a), [mp.mpf(0), lo])
    return rec(0, z)
