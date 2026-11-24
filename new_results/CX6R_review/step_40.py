code = r'''import mpmath as mp, time
import cx6r_common as C
print("CX6R step_13: (15) and (17) actual vs bound at k=99 (identity route + direct pieces)")
mp.mp.dps = 100
tab = C.load_table()
a = C.aa
x0 = 2*mp.exp(mp.pi/4)
K0 = 64/(mp.pi*a)*mp.exp(a/2+mp.pi/4)
T0 = 32/(3*mp.pi*a)*mp.exp(3*a/2+3*mp.pi/4)
def M_main_ray(k, level=6):
    return 2/mp.pi*mp.im(C.glint(lambda r: C.H(r*C.w4)*C.E(r*C.w4)*C.Qk(k,r*C.w4)*C.w4, C.RAY_SEGS, level))
def R2_piece(k, level=6):   # (2/pi) Im int_Gamma H (K-E) Q_k
    def f(r):
        s = mp.mpf('0.5')+r*C.w4
        return C.H(s)*(C.K(s)-C.E(s))*C.Qk(k,s)*C.w4
    return 2/mp.pi*mp.im(C.glint(f, C.GAMMA_SEGS, level))
def R3_piece(k, level=6):   # (2/pi) Im int_0^{1/2} H E Q_k
    return 2/mp.pi*mp.im(C.glint(lambda s: C.H(s)*C.E(s)*C.Qk(k,s), [[0,0.5]], level))
def RHS15(k):
    tot = mp.exp(mp.pi/4)*mp.exp(-mp.sqrt((3*mp.pi/2)*k))
    tot += sum((x0**i/mp.factorial(i))*mp.exp(-mp.sqrt((3*mp.pi/4)*(i+1)*(k-i))) for i in range(1,k+1))
    return K0*tot
k = 99; N = k+1
t=time.time()
mm = M_main_ray(k); R2 = R2_piece(k); R3 = R3_piece(k)
rest = tab[k] - mm
# identity (8): rest = -J_k + R1 + R2 - R3  =>  R1 = rest + J_k - R2 + R3
R1 = rest + C.Jk(k) - R2 + R3
print(f"k={k}: |R1| = |(2/pi) Im int HKA_k| = {mp.nstr(abs(R1),10)}  ({time.time()-t:.0f}s)")
print(f"   RHS(15) = {mp.nstr(RHS15(k),10)}   ratio = {mp.nstr(abs(R1)/RHS15(k),6)}")
print(f"   |R2| = {mp.nstr(abs(R2),10)}   RHS(17) = {mp.nstr(T0*mp.exp(-mp.sqrt(9*mp.pi/4*k)),10)}  ratio = {mp.nstr(abs(R2)/(T0*mp.exp(-mp.sqrt(9*mp.pi/4*k))),6)}")
print(f"   |R3| = {mp.nstr(abs(R3),10)}   (18) bound = {mp.nstr(2/(mp.pi*k)*mp.mpf(3)**(-k),10)}  ratio = {mp.nstr(abs(R3)/(2/(mp.pi*k)*mp.mpf(3)**(-k)),6)}")
print(f"   |J_k| = {mp.nstr(abs(C.Jk(k)),10)}   e^-2sqrtN = {mp.nstr(mp.exp(-2*mp.sqrt(N)),10)}")
print(f"   |rest| = {mp.nstr(abs(rest),10)}  C_*e^-2sqrtN = {mp.nstr((2+T0+K0*(mp.exp(mp.pi/4)+mp.exp(x0)+mp.exp(mp.e*x0)))*mp.exp(-2*mp.sqrt(N)),10)}")
'''
open('step_13.py','w').write(code)
exec(compile(code,'step_13.py','exec'))
