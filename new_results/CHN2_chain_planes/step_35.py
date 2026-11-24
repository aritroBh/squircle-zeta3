
import mpmath as mp, time, pickle
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

u2_ = (mp.mpf('0.3'), mp.mpf('0.8'), mp.mpf('0.5'), mp.mpf('0.6'), mp.mpf('0.4'))
t=time.time(); vt2 = vol_t(*u2_); print("pt2 vol_t:", mp.nstr(vt2,50), f"{time.time()-t:.0f}s", flush=True)
pickle.dump({'vt2':str(vt2)}, open('numcheck_pt2.pkl','wb'))
t=time.time(); vx2 = vol_x(*u2_); print("pt2 vol_x:", mp.nstr(vx2,50), f"{time.time()-t:.0f}s", flush=True)
print("pt2 diff t-x:", mp.nstr(vt2-vx2,5))
pickle.dump({'vt2':str(vt2),'vx2':str(vx2),'diff':str(vt2-vx2)}, open('numcheck_pt2.pkl','wb'))
