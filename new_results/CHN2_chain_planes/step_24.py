
import mpmath as mp
mp.mp.dps = 40

def vol_triple(u1,u2,u3,u4,u5):
    def inner_t4(t3,t2):
        f = lambda t4: u4*t4**(u4-1)*(1-t4)**u5
        return mp.quad(f, [0, 1-t3])
    def inner_t2(t3):
        f = lambda t2: u2*t2**(u2-1)*(1-t2)**u1 * inner_t4(t3,t2)
        return mp.quad(f, [0, 1-t3])
    return mp.quad(lambda t3: u3*t3**(u3-1)*inner_t2(t3), [0,1])

def vol_x(u1,u2,u3,u4,u5):
    p2,p3,p4 = 1/u2,1/u3,1/u4
    def inner_x4(x3):
        lim = (1-x3**p3)**u4
        return mp.quad(lambda x4: (1-x4**p4)**u5, [0, lim])
    def inner_x2(x3):
        lim = (1-x3**p3)**u2
        return mp.quad(lambda x2: (1-x2**p2)**u1 * inner_x4(x3), [0, lim])
    return mp.quad(lambda x3: inner_x2(x3), [0,1])

def vol_dbl(u1,u2,u3,u4,u5, N=120):
    s = mp.mpf(0)
    for n in range(N):
        An = mp.rf(-u1,n)/mp.factorial(n) * u2/(u2+n)
        for m in range(N):
            Bm = mp.rf(-u5,m)/mp.factorial(m) * u4/(u4+m)
            K = mp.gamma(1+u3)*mp.gamma(1+u2+u4+n+m)/mp.gamma(1+u2+u3+u4+n+m)
            s += An*Bm*K
    return s

def vol_A(u1,u2,u3,u4,u5, N=300):
    def A(a1,a2,a4,a5):
        s = mp.mpf(0)
        for m in range(N):
            s += mp.rf(-a5,m)/mp.factorial(m) * a4/(a4+m) * a2*mp.beta(a2+a4+m, 1+a1+u3)
        return s
    return A(u1,u2,u4,u5) + A(u5,u4,u2,u1)

u = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
import time
t=time.time(); v1=vol_triple(*u); print("triple:", mp.nstr(v1, 30), f"{time.time()-t:.1f}s")
t=time.time(); v2=vol_x(*u); print("x-coord:", mp.nstr(v2, 30), f"{time.time()-t:.1f}s")
print("diff:", mp.nstr(v1-v2, 5))
t=time.time(); print("dbl series (120^2, approx):", mp.nstr(vol_dbl(*u), 25), f"{time.time()-t:.1f}s")
t=time.time(); print("A+A_rev (N=300, approx):", mp.nstr(vol_A(*u), 25), f"{time.time()-t:.1f}s")
