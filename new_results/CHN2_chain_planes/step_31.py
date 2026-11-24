
import mpmath as mp, time, pickle
mp.mp.dps = 70
EPS = mp.mpf(10)**(-64)
def ts1d(f, h=mp.mpf('0.07')):
    s = mp.mpf(0)
    for k in range(0, 600):
        t = k*h
        u_ = mp.mpf(mp.pi)/2*mp.sinh(t)
        y = mp.tanh(u_)
        if y >= 1: break
        w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w
        s += term
        if k>25 and abs(term) < EPS: break
    for k in range(1, 600):
        t = -k*h
        u_ = mp.mpf(mp.pi)/2*mp.sinh(t)
        y = mp.tanh(u_)
        if y <= -1: break
        w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w
        s += term
        if k>25 and abs(term) < EPS: break
    return s*h
def int_ab(f,a,b): return ts1d(lambda y: f(a+(y+1)/2*(b-a)))*(b-a)/2

def I2s(T,a1,a2):  # int_0^T a2 t^(a2-1)(1-t)^a1 dt = int_0^{T^a2} (1-w^(1/a2))^a1 dw
    return int_ab(lambda w: (1-w**(1/a2))**a1, 0, T**a2)

def vol_t(u1,u2,u3,u4,u5):
    # outer t3 = s^(1/u3): u3 t3^(u3-1) dt3 = ds
    f = lambda s: I2s(1-s**(1/u3),u1,u2)*I2s(1-s**(1/u3),u5,u4)
    return int_ab(f, 0, mp.mpf('0.5')) + int_ab(f, mp.mpf('0.5'), 1)

u = (mp.mpf('0.5'), mp.mpf(1)/3, mp.mpf('0.25'), mp.mpf('0.2'), mp.mpf(1)/6)
t=time.time(); vt = vol_t(*u); print("vol_t:", mp.nstr(vt,45), f"{time.time()-t:.0f}s")
print("target (vol_x): 0.727963524001991507020626964944895956295200614")
print("diff:", mp.nstr(vt-mp.mpf('0.727963524001991507020626964944895956295200614'),5))
