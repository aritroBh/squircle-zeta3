
import mpmath as mp
mp.mp.dps = 60
u = mp.mpf('0.002')
def pieces(u1,u2,u3,u4,u5, M=20000):
    # returns m0 value, sum_{m>=1} alpha_m for A(u1,u2;u4,u5)
    m0 = u2*mp.beta(u2+u4, 1+u1+u3)
    s = mp.mpf(0)
    t = (-u5)*u4/(u4+1)*u2*mp.beta(u2+u4+1, 1+u1+u3)
    s = t
    S4 = u1+u2+u3+u4
    for m in range(1, M):
        t = t * ((m-u5)/(m+1)) * ((u4+m)/(u4+m+1)) * ((u2+u4+m)/(1+S4+m))
        s += t
    return m0, s
m0A, sA = pieces(0,u,u,u,u)
m0R, sR = pieces(u,u,u,u,0)
print("A   m0 =", mp.nstr(m0A,20), " m>=1 sum =", mp.nstr(sA,20), " /u^3 =", mp.nstr(sA/u**3,12))
print("Arev m0 =", mp.nstr(m0R,20), " m>=1 sum =", mp.nstr(sR,20), " /u^3 =", mp.nstr(sR/u**3,12))
# expected: m>=1 sums ~ -zeta(3)*u^3 each
print("-zeta(3) =", mp.nstr(-mp.zeta(3),12))
# m0 [u^3] extraction
c2A = (m0A - mp.mpf('0.5'))  # m0A = (1/2)(1 + c2 u^2 + c3 u^3...) at a2/(a2+a4)=1/2
print("m0A*u? :", mp.nstr(m0A, 30))
