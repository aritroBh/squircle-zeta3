import mpmath as mp
mp.mp.dps = 25
c0 = mp.mpf('0.31')
B = 2/(1+c0)
def G3(w):
    # G(w1,w2,w3;1) nested quad, letters >1 or 0
    w = [mp.mpf(x) for x in w]
    return mp.quad(lambda t1: mp.quad(lambda t2: mp.quad(lambda t3: 1/(t3-w[2]), [mp.mpf(0),t2])/(t2-w[1]), [mp.mpf(0),t1])/(t1-w[0]), [mp.mpf(0),mp.mpf(1)])
G20B = G3([2,0,B])
G02B = G3([0,2,B])
G0B2 = G3([0,B,2])
GB02 = G3([B,0,2])
G2 = mp.log(1-mp.mpf(1)/2)
GB0 = mp.polylog(2,1/B)   # = G(B,0;1)
lhs = G2*GB0
rhs = -G20B - G02B - G0B2
print("lhs G(2)G(B,0) =", mp.nstr(lhs,20))
print("rhs            =", mp.nstr(rhs,20))
print("diff =", mp.nstr(abs(lhs-rhs),5))
# also un-eliminated: G(2,B,0)+G(B,2,0)+G(B,0,2)
def G_trailing0(u):
    us = [mp.mpf(t) for t in u]
    def rec(idx, lo):
        if idx == len(us)-1:
            return mp.quad(lambda t: mp.log(t)/(t-us[idx]), [mp.mpf(0), lo])
        return mp.quad(lambda t: rec(idx+1, t)/(t-us[idx]), [mp.mpf(0), lo])
    return rec(0, mp.mpf(1))
G2B0 = G_trailing0([2,B]); GB20 = G_trailing0([B,2])
print("direct: G(2,B,0)+G(B,2,0)+G(B,0,2) =", mp.nstr(G2B0+GB20+GB02,20))
print("lhs                                =", mp.nstr(lhs,20))
