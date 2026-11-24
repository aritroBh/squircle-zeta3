code = r'''import mpmath as mp, time
import cx6r_common as C
mp.mp.dps = 60
print("CX6R step_08b: rotation check with nsum extrapolation of period integrals")
def M_main_real_nsum(k):
    cache = {}
    def term(j):
        j = int(j)
        if j not in cache:
            cache[j] = mp.im(C.glint(lambda s: C.H(s)*C.E(s)*C.Qk(k,s), [[4*j, 4*j+4]], 5))
        return cache[j]
    S = mp.nsum(term, [0, mp.inf], tol=mp.mpf('1e-50'))
    return 2/mp.pi*S
def M_main_ray(k, level=6):
    def f(r):
        s = r*C.w4
        return C.H(s)*C.E(s)*C.Qk(k,s)*C.w4
    return 2/mp.pi*mp.im(C.glint(f, C.RAY_SEGS, level))
for k in [3, 10, 30]:
    t=time.time()
    mr = M_main_real_nsum(k)
    my = M_main_ray(k)
    print(f"k={k:3d}  ({time.time()-t:.0f}s)")
    print(f"   real-axis(nsum) = {mp.nstr(mr,30)}")
    print(f"   ray             = {mp.nstr(my,30)}")
    print(f"   |diff| = {mp.nstr(abs(mr-my),5)}")
'''
open('step_08b.py','w').write(code)
exec(compile(code,'step_08b.py','exec'))
