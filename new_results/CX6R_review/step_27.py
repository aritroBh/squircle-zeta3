code = r'''import mpmath as mp, time
import cx6r_common as C
mp.mp.dps = 80
tab = C.load_table()
print("CX6R step_05: eq (7) M_k = -J_k + (2/pi) Im int_Gamma H K P_k ds  vs table")
print("outer: Gamma GL level 6 (96 nodes/seg); inner u-integral GL level 5 (48 nodes/seg)")

def gamma_int_P(k, level):
    def f(r):
        s = mp.mpf('0.5')+r*C.w4
        return C.H(s)*C.K(s)*C.Pk(k, s)*C.w4
    return C.glint(f, C.GAMMA_SEGS, level)

# convergence self-check at k=3
t=time.time()
I5 = gamma_int_P(3, 5)
print(f"k=3 level5 outer: {mp.nstr(I5,30)}  ({time.time()-t:.0f}s)")
t=time.time()
I6 = gamma_int_P(3, 6)
print(f"k=3 level6 outer: {mp.nstr(I6,30)}  ({time.time()-t:.0f}s)")
print("level diff:", mp.nstr(abs(I6-I5),5))
for k in [0,1,2,3,4,5,6,7,8,20]:
    t=time.time()
    I = gamma_int_P(k, 6)
    val = -C.Jk(k) + 2/mp.pi*mp.im(I)
    err = abs(val-tab[k])
    print(f"k={k:3d}  table={mp.nstr(tab[k],24)}")
    print(f"       eq(7)={mp.nstr(val,24)}   |err|={mp.nstr(err,5)}  ({time.time()-t:.0f}s)")
'''
open('step_05.py','w').write(code)
exec(compile(code,'step_05.py','exec'))
