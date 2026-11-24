import mpmath as mp
mp.mp.dps = 25
c0 = mp.mpf('0.31')
B = 2/(1+c0)

def G_trailing0(u):
    us = [mp.mpf(t) for t in u]
    def rec(idx, lo):
        if idx == len(us)-1:
            u = us[idx]
            return mp.quad(lambda t: mp.log(t)/(t-u), [mp.mpf(0), lo])
        u = us[idx]
        return mp.quad(lambda t: rec(idx+1, t)/(t-u), [mp.mpf(0), lo])
    return rec(0, mp.mpf(1))

G2 = mp.log(1-mp.mpf(1)/2)
GB0 = G_trailing0([B])
lhs = G2*GB0
G2B0 = G_trailing0([2,B])
GB20 = G_trailing0([B,2])
# G(B,0,2;1): int dt1/(t1-B) int dt2/t2 int dt3/(t3-2)
GB02 = mp.quad(lambda t1: mp.quad(lambda t2: mp.log(1-t2/2)/t2, [mp.mpf(0),t1])/(t1-B), [mp.mpf(0),mp.mpf(1)])
rhs = G2B0+GB20+GB02
print("G(2) =", mp.nstr(G2,20))
print("G(B,0;1) =", mp.nstr(GB0,20), " Li2(1/B) =", mp.nstr(mp.polylog(2,1/B),20))
print("lhs =", mp.nstr(lhs,20))
print("rhs =", mp.nstr(rhs,20))
print("diff =", mp.nstr(abs(lhs-rhs),5))
