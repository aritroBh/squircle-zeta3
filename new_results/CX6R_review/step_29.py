code = r'''import mpmath as mp, time
import cx6r_common as C
mp.mp.dps = 80
tab = C.load_table()
print("CX6R step_07: eq (8) decomposition + M_main rotation (real axis vs ray)")
def M_main_ray(k, level=6):
    def f(r):
        s = r*C.w4
        return C.H(s)*C.E(s)*C.Qk(k,s)*C.w4
    return 2/mp.pi*mp.im(C.glint(f, C.RAY_SEGS, level))
def M_main_real(k):
    # (2/pi) Im int_0^inf H E Q_k ds on the positive real axis (oscillatory)
    def f(s):
        return C.H(s)*C.E(s)*C.Qk(k,s)
    pts = [0] + [4*j for j in range(1,13)] + [mp.inf]
    return 2/mp.pi*mp.im(mp.quad(f, pts))
def pieces(k, level=6):
    R1 = 2/mp.pi*mp.im(C.glint(lambda r: C.H(mp.mpf('0.5')+r*C.w4)*C.K(mp.mpf('0.5')+r*C.w4)*C.Ak(k, mp.mpf('0.5')+r*C.w4)*C.w4, C.GAMMA_SEGS, level))
    R2 = 2/mp.pi*mp.im(C.glint(lambda r: C.H(mp.mpf('0.5')+r*C.w4)*(C.K(mp.mpf('0.5')+r*C.w4)-C.E(mp.mpf('0.5')+r*C.w4))*C.Qk(k, mp.mpf('0.5')+r*C.w4)*C.w4, C.GAMMA_SEGS, level))
    R3 = 2/mp.pi*mp.im(C.glint(lambda s: C.H(s)*C.E(s)*C.Qk(k,s), [[0,0.5]], level))
    return R1, R2, R3
for k in [3, 10, 30]:
    t=time.time()
    mm_ray = M_main_ray(k)
    mm_real = M_main_real(k)
    R1,R2,R3 = pieces(k)
    rest = -C.Jk(k) + R1 + R2 - R3
    lhs = tab[k] - mm_ray
    print(f"k={k:3d}  ({time.time()-t:.0f}s)")
    print(f"   M_main(ray)  = {mp.nstr(mm_ray,24)}")
    print(f"   M_main(real) = {mp.nstr(mm_real,24)}   |rot err|={mp.nstr(abs(mm_ray-mm_real),5)}")
    print(f"   table - M_main = {mp.nstr(lhs,24)}")
    print(f"   eq(8) rest     = {mp.nstr(rest,24)}   |err|={mp.nstr(abs(rest-lhs),5)}")
    print(f"   pieces: R1={mp.nstr(R1,10)} R2={mp.nstr(R2,10)} R3={mp.nstr(R3,10)} J_k={mp.nstr(C.Jk(k),10)}")
'''
open('step_07.py','w').write(code)
exec(compile(code,'step_07.py','exec'))
