import mpmath as mp
mp.mp.dps = 25
c0 = mp.mpf('0.31')
A = 1/(1+c0); B = 2/(1+c0)

def G_trailing0(u, z=mp.mpf(1)):
    """G(u1..u_{n-1}, 0; z) = int_0^z dt1/(t1-u1) ... int dt_{n-1}/(t_{n-1}-u_{n-1}) * log(t_{n-1})"""
    us = [mp.mpf(t) for t in u]
    def rec(idx, lo):
        if idx == len(us)-1:
            u = us[idx]
            return mp.quad(lambda t: mp.log(t)/(t-u), [mp.mpf(0), lo])
        u = us[idx]
        return mp.quad(lambda t: rec(idx+1, t)/(t-u), [mp.mpf(0), lo])
    return rec(0, mp.mpf(z))

GB0 = G_trailing0([B])
GA0 = G_trailing0([A])
# check G(B,0;1) = -G(0,B;1) = Li_2(1/B)
def Li2(z):
    return mp.polylog(2, z)
print("G(B,0;1) =", mp.nstr(GB0,20), "  Li2(1/B) =", mp.nstr(Li2(1/B),20))
print("G(A,0;1) =", mp.nstr(GA0,20), "  Li2(1/A) =", mp.nstr(Li2(1/A),20))
# equation: G(A;1)*G(B,0;1) =? G(A,B,0)+G(B,A,0)+G(B,0,A)
GA = mp.log(1-1/A); GB = mp.log(1-1/B)
GAB0 = G_trailing0([A,B])
GBA0 = G_trailing0([B,A])
# G(B,0,A;1): good word, no zeros: nested quad
GB0A = mp.quad(lambda t1: mp.quad(lambda t2: mp.log(1-t2/A)/(t2)/(1), [0,t1])/(t1-B), [0,1])
lhs = GA*GB0
rhs = GAB0 + GBA0 + GB0A
print("lhs =", mp.nstr(lhs,20))
print("rhs =", mp.nstr(rhs,20))
print("diff =", mp.nstr(abs(lhs-rhs),5))
