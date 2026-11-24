code = r'''import mpmath as mp, time
import cx6r_common as C
mp.mp.dps = 70
print("CX6R step_08c: wedge deformation int_Gamma HEQ_k = int_ray HEQ_k - int_0^{1/2} HEQ_k")
print("and real-axis = int_Gamma (nsum cross-check)")
def int_Gamma_Q(k):
    return C.glint(lambda r: C.H(mp.mpf('0.5')+r*C.w4)*C.E(mp.mpf('0.5')+r*C.w4)*C.Qk(k,mp.mpf('0.5')+r*C.w4)*C.w4, C.GAMMA_SEGS, 6)
def int_ray_Q(k):
    return C.glint(lambda r: C.H(r*C.w4)*C.E(r*C.w4)*C.Qk(k,r*C.w4)*C.w4, C.RAY_SEGS, 6)
def int_seg_Q(k):
    return C.glint(lambda s: C.H(s)*C.E(s)*C.Qk(k,s), [[0,0.5]], 6)
def M_main_real_nsum(k):
    cache = {}
    def term(j):
        j=int(j)
        if j not in cache:
            cache[j] = mp.im(C.glint(lambda s: C.H(s)*C.E(s)*C.Qk(k,s), [[4*j,4*j+4]], 5))
        return cache[j]
    return 2/mp.pi*mp.nsum(term, [0, mp.inf], tol=mp.mpf('1e-55'))
for k in [3, 10, 30]:
    t=time.time()
    G = int_Gamma_Q(k); R = int_ray_Q(k); S = int_seg_Q(k)
    d1 = abs(G-(R-S))
    print(f"k={k:3d}  |Gamma-(ray-seg)| = {mp.nstr(d1,5)}   ({time.time()-t:.0f}s)")
    mr = M_main_real_nsum(k)
    print(f"      |real-axis - (2/pi)Im Gamma| = {mp.nstr(abs(mr-2/mp.pi*mp.im(G)),5)}")
'''
open('step_08c.py','w').write(code)
exec(compile(code,'step_08c.py','exec'))
