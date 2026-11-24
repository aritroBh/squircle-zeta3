import mpmath as mp
import cx6r_common as C
mp.mp.dps = 90
print("CX6R step_09b: refined t-form quadrature check")
def Ak_t(k, s, level):
    def f(t):
        if t == 0: return mp.mpf(0)
        return t**s*((s*mp.log(t)+mp.log1p(t))**k/(1+t) - (s*mp.log(t))**k)
    pts = [0, mp.mpf('1e-6'), mp.mpf('1e-4'), mp.mpf('1e-3'), mp.mpf('0.02'), mp.mpf('0.2'), mp.mpf('0.6'), mp.mpf(1)]
    segs = [[a,b] for a,b in zip(pts, pts[1:])]
    return (-1)**k/mp.factorial(k)*C.glint(f, segs, level)
s = mp.mpf('0.5')
for k in [0, 5, 10]:
    v_hi = Ak_t(k, s, 7)   # 192 nodes/seg
    v_L  = C.Ak(k, s, level=6)
    print(f"k={k:2d} s=1/2: |t-form(refined) - L-form| = {mp.nstr(abs(v_hi-v_L),5)}")
# and r=1.3, k=10
s = mp.mpf('0.5')+mp.mpf('1.3')*C.w4
v_hi = Ak_t(10, s, 7); v_L = C.Ak(10, s, level=6)
print("k=10 r=1.3 refined:", mp.nstr(abs(v_hi-v_L),5))
