
import mpmath as mp
mp.mp.dps = 60
def A_series(u1,u2,u3,u4,u5, M=4000):
    # A(u1,u2;u4,u5) + A(u5,u4;u2,u1) with Richardson tail
    def one(a1,a2,a4,a5):
        S4 = a1+u3+a5+3
        def alpha0(): return a2*mp.beta(a2+a4, 1+a1+u3)
        def term(m):
            return mp.rf(-a5,m)/mp.factorial(m)*a4/(a4+m)*a2*mp.beta(a2+a4+m, 1+a1+u3)
        s = alpha0() + mp.fsum(term(m) for m in range(1, 2000))
        # tail approx: sum_{m>=2000} alpha_m ~ C m^(1-s)/(s-1)
        return s
    return one(u1,u2,u4,u5)+one(u5,u4,u2,u1)

def single_series(u2,u3,u4,u5, M=4000):
    s = mp.gamma(1+u3)*mp.gamma(1+u2+u4)/mp.gamma(1+u2+u3+u4)
    B = (-u5)*u4/(u4+1)
    for m in range(1,M):
        if m>1: B = B*(m-1-u5)/m*(u4+m-1)/(u4+m)
        K = mp.gamma(1+u3)*mp.gamma(1+u2+u4+m)/mp.gamma(1+u2+u3+u4+m)
        s += B*K
    return s

u = mp.mpf('0.01')
vA = A_series(0,u,u,u,u)
vS = single_series(u,u,u,u)
print("A+A_rev at u1=0:", mp.nstr(vA, 30))
print("single series    :", mp.nstr(vS, 30))
print("diff:", mp.nstr(vA-vS, 10))
c2 = -3*mp.zeta(2); 
print("A: c3 approx:", mp.nstr((vA-1-c2*u**2)/u**3, 12))
print("S: c3 approx:", mp.nstr((vS-1-c2*u**2)/u**3, 12))
