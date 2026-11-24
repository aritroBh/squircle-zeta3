
import mpmath as mp, time
mp.mp.dps = 80
EPS = mp.mpf(10)**(-72)
def make_ts(h):
    h = mp.mpf(h)
    def ts1d(f):
        s = mp.mpf(0)
        for k in range(0, 800):
            t = k*h
            u_ = mp.mpf(mp.pi)/2*mp.sinh(t)
            y = mp.tanh(u_)
            if y >= 1: break
            w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
            term = f(y)*w
            s += term
            if k>30 and abs(term) < EPS: break
        for k in range(1, 800):
            t = -k*h
            u_ = mp.mpf(mp.pi)/2*mp.sinh(t)
            y = mp.tanh(u_)
            if y <= -1: break
            w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
            term = f(y)*w
            s += term
            if k>30 and abs(term) < EPS: break
        return s*h
    return ts1d

def vol_t(u1,u2,u3,u4,u5, h_in='0.05', h_out='0.06'):
    tsi = make_ts(h_in); tso = make_ts(h_out)
    def int_ab_i(f,a,b): return tsi(lambda y: f(a+(y+1)/2*(b-a)))*(b-a)/2
    def int_ab_o(f,a,b): return tso(lambda y: f(a+(y+1)/2*(b-a)))*(b-a)/2
    def I2s(T,a1,a2):
        lim = T**a2
        return int_ab_i(lambda w: (1-w**(1/a2))**a1, 0, lim/2) + int_ab_i(lambda w: (1-w**(1/a2))**a1, lim/2, lim)
    def f(s):
        t3 = s**(1/u3)
        return I2s(1-t3,u1,u2)*I2s(1-t3,u5,u4)
    tot = mp.mpf(0)
    for a,b in [(0,mp.mpf('0.25')),(mp.mpf('0.25'),mp.mpf('0.5')),(mp.mpf('0.5'),mp.mpf('0.75')),(mp.mpf('0.75'),1)]:
        tot += int_ab_o(f, a, b)
    return tot

u = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
t=time.time(); vt = vol_t(*u); print("vol_t:", mp.nstr(vt,50), f"{time.time()-t:.0f}s")
ref = mp.mpf('0.727963524001991507020626964944895956295200614')
print("diff vs vol_x:", mp.nstr(vt-ref,5))
