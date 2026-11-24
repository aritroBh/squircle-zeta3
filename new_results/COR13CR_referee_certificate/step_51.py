import time, sys
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/HEXC2R_referee_certificate')
sys.path.insert(0, base + '/COR7_classical_corners')
import mpmath as mp
import reflib as R
import reflib13 as R13
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)
DPS = 1150
mp.mp.dps = DPS

def my_li(s, z):
    z = mp.mpf(z)
    eps = mp.mpf(10)**(-DPS-30)
    tot = mp.mpf(0); n = 1; zn = z
    while True:
        t = zn/mp.mpf(n)**s
        tot += t
        if n > 200 and abs(t)*abs(z)/(1-abs(z)) < eps: break
        n += 1; zn *= z
    return tot
def my_li4_neg2():
    y = mp.mpc(mp.mpf(1)/2) + mp.log(2)/(2*mp.pi*mp.j)
    B4 = y**4 - 2*y**3 + y**2 - mp.mpf(1)/30
    return (-(2*mp.pi*mp.j)**4/mp.factorial(4) * B4 - my_li(4, mp.mpf(-1)/2)).real

Li21 = my_li(2,"1/3"); Li2m = my_li(2,"-1/3")
Li31 = my_li(3,"1/3"); Li3m = my_li(3,"-1/2")
Li41 = my_li(4,"1/3"); Li4m3 = my_li(4,"-1/3")
Li423 = my_li(4,"2/3"); Li412 = my_li(4,"1/2"); Li419 = my_li(4,"1/9")
Li4m2 = my_li4_neg2()
log(f"Li values done ({time.monotonic()-T0:.0f}s)")
pi, z2, z3 = mp.pi, mp.zeta(2), mp.zeta(3)
L2, L3 = mp.log(2), mp.log(3)

# closed forms (COR7)
W2cl = -Li21 - Li2m + pi**2/12
W3cl = (-26*z3 - 87*Li31 - 130*Li3m + mp.mpf(21)/4*pi**2*L2 - mp.mpf(43)/12*pi**2*L3
        + 67*Li21*L2 - 45*Li21*L3 + mp.mpf(65)/3*L2**3 - 65*L2**2*L3 + mp.mpf(109)/2*L2*L3**2
        - mp.mpf(37)/3*L3**3)
W4cl = (-mp.mpf(22)/3*Li4m2 - mp.mpf(17)/2*Li41 + mp.mpf(20)/3*Li423 + mp.mpf(110)/9*Li412
        + mp.mpf(25)/48*Li419 - 12*Li31*L2 + mp.mpf(28)/3*Li31*L3 - mp.mpf(16)/3*Li3m*L2
        + mp.mpf(16)/3*Li3m*L3 + mp.mpf(14)/3*z3*L2 - mp.mpf(25)/9*z3*L3 - 4*Li21**2
        + mp.mpf(13)/3*Li21*L2**2 - 6*Li21*L2*L3 + mp.mpf(2)/3*Li21*L3**2
        + mp.mpf(11)/18*pi**2*Li21 - mp.mpf(31)/27*pi**2*L2**2 + mp.mpf(31)/72*pi**2*L3**2
        - mp.mpf(271)/1296*pi**4 + mp.mpf(46)/27*L2**4 - mp.mpf(62)/9*L2**3*L3
        + mp.mpf(61)/6*L2**2*L3**2 - mp.mpf(17)/3*L2*L3**3 + mp.mpf(103)/144*L3**4)
M3_closed = W2cl + W3cl + W4cl

# direct assembly from G1 fits with MY hpl_value at 1/3
def fit_val(name, weight):
    tot = mp.mpf(0)
    for e, w, v in R13.fit_terms(name):
        if len(w) != weight: continue
        hv = mp.mpf(1) if not w else R.hpl_value(w, mp.mpf(1)/3)
        tot += mp.mpf(v.numerator)/v.denominator * (mp.mpf(1)/3)**e * hv
    return tot
R3 = sum(fit_val("hpl_fit_k3_(0,0,0).json", w) for w in (2,3,4))
Psi = sum(fit_val("hpl_fit_k3_(1,0,0).json", w) for w in (0,1,2))
M3_hpl = R3 + z2*Psi - z3*(mp.mpf(1)/2 - mp.mpf(2)/3*L2)
log(f"fits assembled ({time.monotonic()-T0:.0f}s)")

# also W-components vs stored references
C7 = base + '/COR7_classical_corners/'
def ref(name): return mp.mpf(open(C7 + name).read().strip())
for name, v, f in [("W2", W2cl, "c13_W2_hpl_1150.txt"), ("W3", W3cl, "c13_W3_hpl_1150.txt"),
                   ("W4", W4cl, "c13_W4_hpl_1150.txt")]:
    d = abs(v - ref(f))
    log(f"{name}: |closed form - {f}| = {mp.nstr(d,3)}")
    assert d < mp.mpf('1e-1080')

quad = open(C7 + 'c13_M3_quad_1100.txt').read().strip()
hpl  = open(C7 + 'c13_M3_hpl_1150.txt').read().strip()
d1 = abs(M3_closed - mp.mpf(quad))
d2 = abs(M3_closed - mp.mpf(hpl))
d3 = abs(M3_hpl - mp.mpf(hpl))
d4 = abs(M3_hpl - M3_closed)
log(f"|M3(closed form) - c13_M3_quad_1100| = {mp.nstr(d1,3)}  ({int(-mp.log10(d1))} digits)")
log(f"|M3(closed form) - c13_M3_hpl_1150| = {mp.nstr(d2,3)}  ({int(-mp.log10(d2))} digits)")
log(f"|M3(my HPL assembly) - c13_M3_hpl_1150| = {mp.nstr(d3,3)} ({int(-mp.log10(d3))} digits)")
log(f"|M3(closed) - M3(my HPL assembly)| = {mp.nstr(d4,3)} ({int(-mp.log10(d4))} digits)")
assert d1 < mp.mpf('1e-1080') and d2 < mp.mpf('1e-1080') and d3 < mp.mpf('1e-1080') and d4 < mp.mpf('1e-1080')
log(f"M3 = {mp.nstr(M3_closed, 60)}")
with open('progress.log','a') as f:
    f.write("referee step(M3 1150dps): closed form vs quad(1100) %d digits, vs hpl(1150) %d; my assembly agrees %d digits\n"
            % (int(-mp.log10(d1)), int(-mp.log10(d2)), int(-mp.log10(d4))))
log("STEP PASS")
