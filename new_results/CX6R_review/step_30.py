code = r'''import mpmath as mp, time
import cx6r_common as C
mp.mp.dps = 60
print("CX6R step_08: rotation check M_main(real axis) = M_main(ray), sec.3")
def period_int(k, j):
    # integral over [4j, 4j+4] (E has period 4 on the real axis)
    return C.glint(lambda s: C.H(s)*C.E(s)*C.Qk(k,s), [[4*j, 4*j+4]], 5)
def M_main_real_careful(k, T=400):
    per = {}
    def term(j):
        if j not in per: per[j] = period_int(k, j)
        return per[j]
    # integrand ~ s^{-5/2} oscillatory; period integrals decrease ~ j^{-5/2}: nsum OK
    S = mp.nsum(term, [0, int(T/4)-1])
    # tail beyond T: |int_T^inf| <= (4/pi) * sqrt(pi/2) T^{-3/2} * (T/(T+1))^k/(T+1)  (parts)
    tail = (4/mp.pi)*mp.sqrt(mp.pi/2)*T**(-mp.mpf(3)/2)*(mp.mpf(T)/(T+1))**k/(T+1)
    return 2/mp.pi*mp.im(S), 2/mp.pi*tail
def M_main_ray(k, level=6):
    def f(r):
        s = r*C.w4
        return C.H(s)*C.E(s)*C.Qk(k,s)*C.w4
    return 2/mp.pi*mp.im(C.glint(f, C.RAY_SEGS, level))
for k in [3, 10, 30]:
    t=time.time()
    mr, tail = M_main_real_careful(k)
    my = M_main_ray(k)
    print(f"k={k:3d}  ({time.time()-t:.0f}s)")
    print(f"   real-axis = {mp.nstr(mr,26)}  (tail allowance {mp.nstr(tail,5)})")
    print(f"   ray       = {mp.nstr(my,26)}")
    print(f"   |diff| = {mp.nstr(abs(mr-my),5)}")
'''
open('step_08.py','w').write(code)
exec(compile(code,'step_08.py','exec'))
