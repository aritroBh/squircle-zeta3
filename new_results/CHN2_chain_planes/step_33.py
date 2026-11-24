
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

u = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
t=time.time(); vx = vol_x(*u); print("vol_x:", mp.nstr(vx,50), f"{time.time()-t:.0f}s")
vt = mp.mpf('0.72796352400199149975029831223646699601899153661546')
print("vx - vt:", mp.nstr(vx-vt,5))
