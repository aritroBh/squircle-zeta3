import mpmath as mp, time
import cx6r_common as C
print("CX6R step_17: more worst-case sweeps + k=100/200/400 rest vs arb log")
mp.mp.dps = 100
tab = C.load_table()
a = C.aa
# (15) first piece F(k,2;a) termwise
def Factual(j, q, lam, level=7):
    rs = mp.sqrt(j*q/(mp.sqrt(2)*lam)) if j > 0 else mp.mpf(1)
    pts = [0, rs/4, rs/2, rs, 2*rs, 4*rs, 10*rs, 10*rs+40/lam]
    segs = [[x,y] for x,y in zip(pts, pts[1:])]
    def f(r):
        s = mp.mpf('0.5')+r*C.w4
        return mp.exp(-lam*r)*abs(s)**j/abs(s+q)**(j+1)
    return C.glint(f, segs, level)
def RHS14(j, q, lam):
    return 4/lam*mp.exp(lam/mp.mpf(2)+lam*q/mp.sqrt(2))*mp.exp(-mp.sqrt(3*lam*j*q/mp.sqrt(2)))
for k in (99, 200, 400):
    print(f"F({k},2,a)/RHS14 = {mp.nstr(Factual(k,2,a,7)/RHS14(k,2,a),6)}")
# (15),(17) actual via identity route at k=200,400
x0 = 2*mp.exp(mp.pi/4)
K0 = 64/(mp.pi*a)*mp.exp(a/2+mp.pi/4)
T0 = 32/(3*mp.pi*a)*mp.exp(3*a/2+3*mp.pi/4)
def M_main_ray(k, level=6):
    return 2/mp.pi*mp.im(C.glint(lambda r: C.H(r*C.w4)*C.E(r*C.w4)*C.Qk(k,r*C.w4)*C.w4, C.RAY_SEGS, level))
def RHS15(k):
    tot = mp.exp(mp.pi/4)*mp.exp(-mp.sqrt((3*mp.pi/2)*k))
    tot += sum((x0**i/mp.factorial(i))*mp.exp(-mp.sqrt((3*mp.pi/4)*(i+1)*(k-i))) for i in range(1,k+1))
    return K0*tot
logvals = {100: '1.69252447066814635e-13', 200: '3.81140793878010622e-18', 400: '4.79912536555363781e-25'}
for k in (100, 200, 400):
    mm = M_main_ray(k)
    rest = abs(tab[k]-mm)
    N = k+1
    Cstar = 2 + T0 + K0*(mp.exp(mp.pi/4)+mp.exp(x0)+mp.exp(mp.e*x0))
    print(f"k={k}: my |M_k-M_main| = {mp.nstr(rest,18)}  (arb log: {logvals[k]})")
    print(f"      rest/(C_*e^-2sqrtN) = {mp.nstr(rest/(Cstar*mp.exp(-2*mp.sqrt(N))),6)}")
# (15) actual at k=200,400 via identity
def R2_piece(k, level=6):
    def f(r):
        s = mp.mpf('0.5')+r*C.w4
        return C.H(s)*(C.K(s)-C.E(s))*C.Qk(k,s)*C.w4
    return 2/mp.pi*mp.im(C.glint(f, C.GAMMA_SEGS, level))
def R3_piece(k, level=6):
    return 2/mp.pi*mp.im(C.glint(lambda s: C.H(s)*C.E(s)*C.Qk(k,s), [[0,0.5]], level))
for k in (200, 400):
    mm = M_main_ray(k); R2 = R2_piece(k); R3 = R3_piece(k)
    rest = tab[k]-mm
    R1 = rest + C.Jk(k) - R2 + R3
    print(f"k={k}: |R1|={mp.nstr(abs(R1),8)} RHS15={mp.nstr(RHS15(k),8)} ratio={mp.nstr(abs(R1)/RHS15(k),6)}")
    print(f"       |R2|={mp.nstr(abs(R2),8)} RHS17={mp.nstr(T0*mp.exp(-mp.sqrt(9*mp.pi/4*k)),8)} ratio={mp.nstr(abs(R2)/(T0*mp.exp(-mp.sqrt(9*mp.pi/4*k))),6)}")
