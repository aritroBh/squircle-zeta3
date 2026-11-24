
import mpmath as mp
mp.mp.dps = 40

def vol_triple(u1,u2,u3,u4,u5):
    # t-integral: u2u3u4 ∫∫∫_{t2+t3<=1, t3+t4<=1} t2^(u2-1)t3^(u3-1)t4^(u4-1)(1-t2)^u1(1-t4)^u5
    f = lambda t3,t2,t4: u2*u3*u4 * t2**(u2-1)*t3**(u3-1)*t4**(u4-1)*(1-t2)**u1*(1-t4)**u5
    return mp.quad(f, [0,1], [0, lambda t3: 1-t3], [0, lambda t3,t2: 1-t3])

def vol_x(u1,u2,u3,u4,u5):
    # direct x-coordinate integral (orthant): dx2 dx3 dx4 with (1-x2^p2)^u1 (1-x4^p4)^u5
    p2,p3,p4 = 1/u2,1/u3,1/u4
    f = lambda x3,x2,x4: (1-x2**p2)**u1*(1-x4**p4)**u5
    return mp.quad(f, [0,1],[0, lambda x3: (1-x3**p3)**u2],[0, lambda x3,x2: (1-x3**p3)**u4])

def vol_dbl(u1,u2,u3,u4,u5, N=80):
    s = mp.mpf(0)
    for n in range(N):
        An = mp.rf(-u1,n)/mp.factorial(n) * u2/(u2+n)
        for m in range(N):
            Bm = mp.rf(-u5,m)/mp.factorial(m) * u4/(u4+m)
            K = mp.gamma(1+u3)*mp.gamma(1+u2+u4+n+m)/mp.gamma(1+u2+u3+u4+n+m)
            s += An*Bm*K
    return s  # tail ~ 1/N^... only approximate; use for rough check

def vol_A(u1,u2,u3,u4,u5, N=200):
    def A(a1,a2,a4,a5):
        s = mp.mpf(0)
        for m in range(N):
            s += mp.rf(-a5,m)/mp.factorial(m) * a4/(a4+m) * a2*mp.beta(a2+a4+m, 1+a1+u3)
        return s
    return A(u1,u2,u4,u5) + A(u5,u4,u2,u1)

u = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
print("triple:", mp.nstr(vol_triple(*u), 30))
print("x-coord:", mp.nstr(vol_x(*u), 30))
print("dbl series (80x80, approx):", mp.nstr(vol_dbl(*u), 25))
print("A+A_rev (N=200, approx):", mp.nstr(vol_A(*u), 25))
