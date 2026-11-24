
import mpmath as mp, time
mp.mp.dps = 80
EPS = mp.mpf(10)**(-72)
def make_ts(h):
    h = mp.mpf(h)
    def ts1d(f):
        s = mp.mpf(0)
        for k in range(0, 800):
            t = k*h; u_ = mp.mpf(mp.pi)/2*mp.sinh(t); y = mp.tanh(u_)
            if y >= 1: break
            w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
            term = f(y)*w; s += term
            if k>30 and abs(term) < EPS: break
        for k in range(1, 800):
            t = -k*h; u_ = mp.mpf(mp.pi)/2*mp.sinh(t); y = mp.tanh(u_)
            if y <= -1: break
            w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
            term = f(y)*w; s += term
            if k>30 and abs(term) < EPS: break
        return s*h
    return ts1d
tsi = make_ts('0.05'); tso = make_ts('0.06')
def int_ab_i(f,a,b): return tsi(lambda y: f(a+(y+1)/2*(b-a)))*(b-a)/2
def int_ab_o(f,a,b): return tso(lambda y: f(a+(y+1)/2*(b-a)))*(b-a)/2
def I2s(T,a1,a2):
    lim = T**a2
    return int_ab_i(lambda w: (1-w**(1/a2))**a1, 0, lim/2) + int_ab_i(lambda w: (1-w**(1/a2))**a1, lim/2, lim)
def vol_t(u1,u2,u3,u4,u5):
    def f(s):
        t3 = s**(1/u3)
        return I2s(1-t3,u1,u2)*I2s(1-t3,u5,u4)
    tot = mp.mpf(0)
    for a,b in [(0,mp.mpf('0.25')),(mp.mpf('0.25'),mp.mpf('0.5')),(mp.mpf('0.5'),mp.mpf('0.75')),(mp.mpf('0.75'),1)]:
        tot += int_ab_o(f, a, b)
    return tot
def vol_x(u1,u2,u3,u4,u5):
    p2,p3,p4 = 1/u2,1/u3,1/u4
    def J2(x3):
        L2 = (1-x3**p3)**u2
        return int_ab_i(lambda x2: (1-x2**p2)**u1, 0, L2/2)+int_ab_i(lambda x2: (1-x2**p2)**u1, L2/2, L2)
    def J4(x3):
        L4 = (1-x3**p3)**u4
        return int_ab_i(lambda x4: (1-x4**p4)**u5, 0, L4/2)+int_ab_i(lambda x4: (1-x4**p4)**u5, L4/2, L4)
    tot = mp.mpf(0)
    for a,b in [(0,mp.mpf('0.25')),(mp.mpf('0.25'),mp.mpf('0.5')),(mp.mpf('0.5'),mp.mpf('0.75')),(mp.mpf('0.75'),1)]:
        tot += int_ab_o(lambda x3: J2(x3)*J4(x3), a, b)
    return tot

def dbl_resum(u1,u2,u3,u4,u5, M=1200, K=6):
    # double series summed over n via 3F2, then over m with tail Richardson
    S4 = u1+u2+u3+u4
    def term(m):
        c = mp.rf(-u5,m)/mp.factorial(m)*u4/(u4+m)
        g = mp.gamma(1+u2+u4+m)/mp.gamma(1+S4+m)
        h = mp.hyper([-u1,u2,1+u2+u4+m],[1+u2,1+S4+m],1)
        return mp.gamma(1+u3)*c*g*h
    Ms = sorted(set([M, M*3//4, M//2, M*3//8, M//4, M//5, M//7]))
    s_exp = u3+u5+2   # term ~ m^-(u3+u5+2)
    vals=[]
    for Mi in Ms:
        vals.append(mp.fsum(term(m) for m in range(Mi)))
    n=len(Ms)
    A = mp.matrix(n, n)
    for i,Mi in enumerate(Ms):
        A[i,0] = 1
        for j in range(1,n):
            A[i,j] = Mi**(1-s_exp) * Mi**(-(j-1))
    return mp.lu_solve(A, mp.matrix(vals))[0]

u2_ = (mp.mpf('0.3'), mp.mpf('0.8'), mp.mpf('0.5'), mp.mpf('0.6'), mp.mpf('0.4'))
t=time.time(); vt2 = vol_t(*u2_); print("pt2 vol_t:", mp.nstr(vt2,45), f"{time.time()-t:.0f}s")
t=time.time(); vx2 = vol_x(*u2_); print("pt2 vol_x:", mp.nstr(vx2,45), f"{time.time()-t:.0f}s")
print("pt2 diff t-x:", mp.nstr(vt2-vx2,5))
t=time.time(); vd = dbl_resum(*u2_); print("pt2 dbl(3F2 resum):", mp.nstr(vd,45), f"{time.time()-t:.0f}s")
print("dbl - vol_t:", mp.nstr(vd-vt2,5))
# anchor for vol_x
print("anchor vol_x(1,1,1,1,1):", mp.nstr(vol_x(1,1,1,1,1),25), "vs 2/15:", mp.nstr(mp.mpf(2)/15,25))
