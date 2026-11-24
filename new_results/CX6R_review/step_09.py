import mpmath as mp, time
import cx6r_common as C
mp.mp.dps = 70
print("CX6R step_09: eq (10) L-form of A_k vs t-form P_k - Q_k")
def Ak_t(k, s):
    def f(t):
        if t == 0: return mp.mpf(0)
        return t**s*((s*mp.log(t)+mp.log1p(t))**k/(1+t) - (s*mp.log(t))**k)
    pts = [0, mp.mpf('0.02'), mp.mpf('0.2'), mp.mpf('0.6'), mp.mpf(1)]
    segs = [[a,b] for a,b in zip(pts, pts[1:])]
    return (-1)**k/mp.factorial(k)*C.glint(f, segs, 6)
for k in [0, 1, 2, 5, 10]:
    for r in [mp.mpf('0'), mp.mpf('1.3'), mp.mpf('7')]:
        s = mp.mpf('0.5')+r*C.w4
        v1 = Ak_t(k, s); v2 = C.Ak(k, s, level=6)
        print(f"k={k:2d} r={mp.nstr(r,4)}: |t-form - L-form| = {mp.nstr(abs(v1-v2),5)}")
