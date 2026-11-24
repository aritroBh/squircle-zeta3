
import mpmath as mp, time, pickle
mp.mp.dps = 50

def ts1d(f, h=mp.mpf('0.09'), tmax=mp.mpf('5.5')):
    s = mp.mpf(0)
    for k in range(0, 400):
        t = k*h
        u_ = mp.mpf(mp.pi)/2*mp.sinh(t)
        y = mp.tanh(u_)
        if y >= 1: break
        w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w
        s += term
        if k>20 and abs(term) < mp.mpf(10)**(-48): break
    for k in range(1, 400):
        t = -k*h
        u_ = mp.mpf(mp.pi)/2*mp.sinh(t)
        y = mp.tanh(u_)
        if y <= -1: break
        w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w
        s += term
        if k>20 and abs(term) < mp.mpf(10)**(-48): break
    return s*h

def int01(f): return ts1d(lambda y: f((y+1)/2))/2

def I2(T,a1,a2): return T**a2 * mp.hyper([a2,-a1],[1+a2], T)

def vol_t(u1,u2,u3,u4,u5):
    f = lambda t3: u3*t3**(u3-1) * I2(1-t3,u1,u2)*I2(1-t3,u5,u4)
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

def alpha_sum(u1,u2,u3,u4,u5, M):
    S4 = u1+u2+u3+u4
    a = (-u5)*u4/(u4+1)*u2*mp.beta(u2+u4+1, 1+u1+u3)
    s = a
    for m in range(1, M):
        a = a * ((m-u5)/(m+1)) * ((u4+m)/(u4+m+1)) * ((u2+u4+m)/(1+S4+m))
        s += a
    return s

def richardson(Ms, vals, s, K):
    n = len(Ms)
    A = mp.matrix(n, K+1)
    for i,M in enumerate(Ms):
        A[i,0] = 1
        for j in range(1,K+1):
            A[i,j] = M**(1-s) * M**(-(j-1))
    b = mp.matrix(vals)
    return mp.lu_solve(A, b)[0]

def vol_A(u1,u2,u3,u4,u5, M=2400, K=5):
    def one(a1,a2,a4,a5):
        s = a1+u3+a5+3
        Ms = [M>>i for i in (K, K-1, K-2, K-1, K-2//2, 0)]
        Ms = sorted(set([m for m in Ms if m>0]))
        vals = [alpha_sum(a1,a2,u3,a4,a5, Mi) for Mi in Ms]
        m0 = a2*mp.beta(a2+a4, 1+a1+u3)
        return m0 + richardson([mp.mpf(x) for x in Ms], vals, s, len(Ms)-1)
    return one(u1,u2,u4,u5) + one(u5,u4,u2,u1)

u = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
t=time.time(); vt = vol_t(*u); print("vol_t:", mp.nstr(vt,40), f"{time.time()-t:.0f}s")
t=time.time(); vx = vol_x(*u); print("vol_x:", mp.nstr(vx,40), f"{time.time()-t:.0f}s")
t=time.time(); va = vol_A(*u); print("vol_A:", mp.nstr(va,40), f"{time.time()-t:.0f}s")
print("vt-vx:", mp.nstr(vt-vx,6), " vt-va:", mp.nstr(vt-va,6))
pickle.dump({'vt':str(vt),'vx':str(vx),'va':str(va)}, open('numcheck_pt1.pkl','wb'))
