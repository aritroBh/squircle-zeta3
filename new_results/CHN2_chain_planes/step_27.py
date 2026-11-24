
import mpmath as mp, time
mp.mp.dps = 50

def ts1d(f, h=mp.mpf('0.09'), tmax=mp.mpf('4.5')):
    # int_{-1}^1 f(y) dy via tanh-sinh
    s = mp.mpf(0); k = 0
    while True:
        t = k*h
        sh = mp.sinh(t); u_ = mp.mpf(mp.pi)/2*sh
        y = mp.tanh(u_); w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w
        s += term
        if k>0 and abs(term) < mp.mpf(10)**(-45) and t>1: break
        k += 1
        if t > tmax and k>10: break
    # symmetric
    s2 = mp.mpf(0); k=1
    while True:
        t = -k*h
        sh = mp.sinh(t); u_ = mp.mpf(mp.pi)/2*sh
        y = mp.tanh(u_); w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w
        s2 += term
        if abs(term) < mp.mpf(10)**(-45) and t<-1: break
        k += 1
        if t < -tmax: break
    return (s+s2)*h

def int01(f):
    return ts1d(lambda y: f((y+1)/2))/2

def I2(T,u1,u2):  # int_0^T u2 t2^(u2-1)(1-t2)^u1 dt2
    return T**u2 * mp.hyper([u2,-u1],[1+u2], T)

def vol_t(u1,u2,u3,u4,u5):
    f = lambda t3: u3*t3**(u3-1) * I2(1-t3,u1,u2)*I2(1-t3,u5,u4)
    # split [0,1] at 1/2 for the endpoint singularity at 0
    v1 = int01(lambda x: f(x/2))/2
    v2 = int01(lambda x: f((x+1)/2))/2
    return v1+v2

def vol_x(u1,u2,u3,u4,u5):
    p2,p3,p4 = 1/u2,1/u3,1/u4
    def J2(x3):
        L2 = (1-x3**p3)**u2
        return int01(lambda y: (1-(y*L2)**p2)**u1)*L2
    def J4(x3):
        L4 = (1-x3**p3)**u4
        return int01(lambda y: (1-(y*L4)**p4)**u5)*L4
    return int01(lambda x3: J2(x3)*J4(x3))

def alpha_series(u1,u2,u3,u4,u5, M):
    S4 = u1+u2+u3+u4
    a = u2*mp.beta(u2+u4, 1+u1+u3)  # m=0
    tot_m0 = a
    # m=1
    a = (-u5)*u4/(u4+1)*u2*mp.beta(u2+u4+1, 1+u1+u3)
    s = a
    arr=[a]
    for m in range(1, M):
        a = a * ((m-u5)/(m+1)) * ((u4+m)/(u4+m+1)) * ((u2+u4+m)/(1+S4+m))
        s += a; arr.append(a)
    return tot_m0, s

def richardson(vals, Ms, s):
    # vals = S(M) at Ms; fit Sinf + M^(1-s)(c0 + c1/M + c2/M^2 + c3/M^3)
    n = len(Ms)
    A = mp.matrix(n, 4)
    for i,M in enumerate(Ms):
        A[i,0] = 1
        for j in range(1,4):
            A[i,j] = M**(1-s) * M**(-(j-1))
    b = mp.matrix(vals)
    sol = mp.lu_solve(A, b)
    return sol[0]

def vol_A(u1,u2,u3,u4,u5, M=1600):
    def one(a1,a2,a4,a5):
        s = u1+u3+u5+3
        Ms = [M//8, M//4, M//2, M]
        vals=[]
        for Mi in Ms:
            _, v = alpha_series(a1,a2,u3,a4,a5, Mi)
            vals.append(v)
        m0,_ = alpha_series(a1,a2,u3,a4,a5, 0)
        return m0 + richardson(vals, [mp.mpf(x) for x in Ms], s)
    return one(u1,u2,u4,u5) + one(u5,u4,u2,u1)

u=(mp.mpf(1),)*5
print("u=(1,1,1,1,1): vol_t =", mp.nstr(vol_t(*u),30), " exact 2/15 =", mp.nstr(mp.mpf(2)/15,30))
t=time.time()
u = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
vt = vol_t(*u); print("vol_t:", mp.nstr(vt,35), f"{time.time()-t:.0f}s")
t=time.time(); vx = vol_x(*u); print("vol_x:", mp.nstr(vx,35), f"{time.time()-t:.0f}s")
print("vol_t - vol_x:", mp.nstr(vt-vx,6))
t=time.time(); va = vol_A(*u); print("vol_A:", mp.nstr(va,35), f"{time.time()-t:.0f}s")
print("vol_t - vol_A:", mp.nstr(vt-va,6))
