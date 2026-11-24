import mpmath as mp
mp.mp.dps = 50

def Li_ab(a, b, x, y, terms=400):
    # Li_{a,b}(x,y) = sum_{0<k1<k2} x^k1 y^k2 / (k1^a k2^b)
    x = mp.mpf(x); y = mp.mpf(y)
    s = mp.mpf(0)
    k2 = 2
    # sum until tail small
    term = mp.mpf(1)
    k1max = terms
    for k2v in range(2, terms):
        inner = mp.mpf(0)
        for k1v in range(1, k2v):
            inner += x**k1v / mp.power(k1v, a)
        s += y**k2v / mp.power(k2v, b) * inner
    return s

def G_int(args, z):
    # G(a1,...,an; z) = int_0^z dt1/(t1-a1) int_0^{t1} dt2/(t2-a2) ...
    args = [mp.mpf(a) for a in args]
    z = mp.mpf(z)
    def rec(idx, lo):
        if idx == len(args):
            return mp.mpf(1)
        a = args[idx]
        f = lambda t: rec(idx+1, t)/(t-a)
        return mp.quad(f, [mp.mpf(0), lo])
    return rec(0, z)

x, y = mp.mpf('0.3'), mp.mpf('0.4')
L11 = Li_ab(1,1,x,y)
G1 = G_int([1/y, 1/(x*y)], 1)   # candidate A: a=1/y, b=1/(xy)
G2 = G_int([1/x, 1/(x*y)], 1)   # candidate B: a=1/x
print("Li_11(x,y) series :", mp.nstr(L11, 30))
print("G(1/y,1/(xy);1)   :", mp.nstr(G1, 30))
print("G(1/x,1/(xy);1)   :", mp.nstr(G2, 30))
print("Li_11(y,x) series :", mp.nstr(Li_ab(1,1,y,x), 30))
