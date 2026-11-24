
import mpmath as mp, time, pickle
mp.mp.dps = 70
EPS = mp.mpf(10)**(-62)

def ts1d(f, h=mp.mpf('0.075')):
    s = mp.mpf(0)
    for k in range(0, 500):
        t = k*h
        u_ = mp.mpf(mp.pi)/2*mp.sinh(t)
        y = mp.tanh(u_)
        if y >= 1: break
        w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w
        s += term
        if k>25 and abs(term) < EPS: break
    for k in range(1, 500):
        t = -k*h
        u_ = mp.mpf(mp.pi)/2*mp.sinh(t)
        y = mp.tanh(u_)
        if y <= -1: break
        w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w
        s += term
        if k>25 and abs(term) < EPS: break
    return s*h

def int01(f): return ts1d(lambda y: f((y+1)/2))/2
def int_ab(f,a,b): return ts1d(lambda y: f(a+(y+1)/2*(b-a)))*(b-a)/2

def I2q(T,a1,a2):  # int_0^T a2 t^(a2-1)(1-t)^a1 dt, pure tanh-sinh
    return int_ab(lambda t: a2*t**(a2-1)*(1-t)**a1, 0, T)

def vol_t(u1,u2,u3,u4,u5):
    def f(t3):
        return u3*t3**(u3-1) * I2q(1-t3,u1,u2)*I2q(1-t3,u5,u4)
    return int_ab(f, 0, mp.mpf('0.5')) + int_ab(f, mp.mpf('0.5'), 1)

def vol_x(u1,u2,u3,u4,u5):
    p2,p3,p4 = 1/u2,1/u3,1/u4
    def J2(x3):
        L2 = (1-x3**p3)**u2
        return int_ab(lambda x2: (1-x2**p2)**u1, 0, L2)
    def J4(x3):
        L4 = (1-x3**p3)**u4
        return int_ab(lambda x4: (1-x4**p4)**u5, 0, L4)
    return int_ab(lambda x3: J2(x3)*J4(x3), 0, 1)

def alpha_sum(u1,u2,u3,u4,u5, M):
    S4 = u1+u2+u3+u4
    a = (-u5)*u4/(u4+1)*u2*mp.beta(u2+u4+1, 1+u1+u3)
    s = a
    for m in range(1, M):
        a = a * ((m-u5)/(m+1)) * ((u4+m)/(u4+m+1)) * ((u2+u4+m)/(1+S4+m))
        s += a
    return s

def vol_A(u1,u2,u3,u4,u5, M=3000):
    def one(a1,a2,a4,a5):
        s = a1+u3+a5+3
        Ms = [M, M*3//4, M//2, M*3//8, M//4, M//6]
        vals = [alpha_sum(a1,a2,u3,a4,a5, Mi) for Mi in Ms]
        K = len(Ms)-1
        A = mp.matrix(len(Ms), K+1)
        for i,Mi in enumerate(Ms):
            A[i,0] = 1
            for j in range(1,K+1):
                A[i,j] = Mi**(1-s) * Mi**(-(j-1))
        b = mp.matrix(vals)
        sol = mp.lu_solve(A, b)
        m0 = a2*mp.beta(a2+a4, 1+a1+u3)
        return m0 + sol[0]
    return one(u1,u2,u4,u5) + one(u5,u4,u2,u1)

u = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
t=time.time(); vt = vol_t(*u); print("vol_t:", mp.nstr(vt,45), f"{time.time()-t:.0f}s")
t=time.time(); vx = vol_x(*u); print("vol_x:", mp.nstr(vx,45), f"{time.time()-t:.0f}s")
t=time.time(); va = vol_A(*u); print("vol_A:", mp.nstr(va,45), f"{time.time()-t:.0f}s")
print("vt-vx:", mp.nstr(vt-vx,5))
print("vt-va:", mp.nstr(vt-va,5))
