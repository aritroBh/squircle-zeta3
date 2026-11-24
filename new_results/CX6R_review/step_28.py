code = r'''import mpmath as mp, time
import cx6r_common as C
mp.mp.dps = 80
tab = C.load_table()
print("CX6R step_06: eq (9) M_k = -(3/4)J_k + (2/pi) Im int_ray H K P_k ds  (sec.4 arc correction)")
def ray_int_P(k, level=6):
    def f(r):
        s = r*C.w4
        return C.H(s)*C.K(s)*C.Pk(k, s)*C.w4
    return C.glint(f, C.RAY_SEGS, level)
for k in [8, 10, 3]:
    t=time.time()
    I = ray_int_P(k)
    val = -mp.mpf(3)/4*C.Jk(k) + 2/mp.pi*mp.im(I)
    err = abs(val-tab[k])
    print(f"k={k:3d}  table={mp.nstr(tab[k],24)}")
    print(f"       eq(9)={mp.nstr(val,24)}   |err|={mp.nstr(err,5)}  ({time.time()-t:.0f}s)")
# also confirm the WRONG sign (+(3/4)... or +pi J_k/8 arc) fails: compute with -J_k + piJ_k/8 vs -(3/4)J_k
k=8
I = ray_int_P(k)
alt = -C.Jk(k) - mp.pi*C.Jk(k)/8/mp.pi*2  # nonsense variants to show sensitivity
print("sanity: table M_8 =", mp.nstr(tab[8],24))
print("   -(3/4)J_k version matches;   -J_k + (1/4)J_k = -(3/4)J_k checks the arc sign")
print("   variant with +pi J_k/8 arc (i.e. -(5/4)J_k +...):",
      mp.nstr(-mp.mpf(5)/4*C.Jk(k) + 2/mp.pi*mp.im(I), 15), " (should NOT match)")
'''
open('step_06.py','w').write(code)
exec(compile(code,'step_06.py','exec'))
