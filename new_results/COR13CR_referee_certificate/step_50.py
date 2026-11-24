import time, sys
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/HEXC2R_referee_certificate')
import mpmath as mp
import reflib as R
import reflib13 as R13
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)
DPS = 250
mp.mp.dps = DPS

# --- my own Li values (direct series; Jonquiere for Li4(-2)) ---
def my_li(s, z):
    z = mp.mpf(z)
    assert abs(z) < 1
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
    rhs = -(2*mp.pi*mp.j)**4/mp.factorial(4) * B4
    return (rhs - my_li(4, mp.mpf(-1)/2)).real
# cross-check vs mpmath
cc = max(max(abs(my_li(s,z) - mp.polylog(s, mp.mpf(z))) for s,z in
             [(2,"1/3"),(2,"-1/3"),(3,"1/3"),(3,"-1/2"),(4,"1/3"),(4,"-1/3"),(4,"2/3"),(4,"1/2"),(4,"-1/2"),(4,"1/9")]),
         abs(my_li4_neg2() - mp.polylog(4, -2)))
log(f"my Li values vs mpmath polylog (10 points + Jonquiere Li4(-2)): worst {mp.nstr(cc,3)}")
assert cc < mp.mpf('1e-230')

# --- my own HPL-based R3^(w)(1/3) and Psi^(w)(1/3) from the G1 fits ---
def fit_val(name, weight):
    tot = mp.mpf(0)
    for e, w, v in R13.fit_terms(name):
        if len(w) != weight: continue
        hv = mp.mpf(1) if not w else R.hpl_value(w, mp.mpf(1)/3)
        tot += mp.mpf(v.numerator)/v.denominator * (mp.mpf(1)/3)**e * hv
    return tot
R32 = fit_val("hpl_fit_k3_(0,0,0).json", 2)
R33 = fit_val("hpl_fit_k3_(0,0,0).json", 3)
R34 = fit_val("hpl_fit_k3_(0,0,0).json", 4)
P0  = fit_val("hpl_fit_k3_(1,0,0).json", 0)
P1  = fit_val("hpl_fit_k3_(1,0,0).json", 1)
P2  = fit_val("hpl_fit_k3_(1,0,0).json", 2)
log(f"R3^(2,3,4)(1/3) and Psi^(0,1,2)(1/3) evaluated with MY hpl_value")
# Psi^(0) must be 1/2; Psi^(1) = (21/2) log 2 - (15/2) log 3 (COR7)
log(f"|Psi^(0) - 1/2| = {mp.nstr(abs(P0 - mp.mpf(1)/2),3)}")
log(f"|Psi^(1) - ((21/2)log2 - (15/2)log3)| = {mp.nstr(abs(P1 - (mp.mpf(21)/2*mp.log(2) - mp.mpf(15)/2*mp.log(3))),3)}")

pi = mp.pi; z2 = mp.zeta(2); z3 = mp.zeta(3)
L2, L3 = mp.log(2), mp.log(3)
Li21, Li2m = my_li(2,"1/3"), my_li(2,"-1/3")
Li31, Li3m = my_li(3,"1/3"), my_li(3,"-1/2")
Li41, Li4m3, Li423, Li412, Li419 = my_li(4,"1/3"), my_li(4,"-1/3"), my_li(4,"2/3"), my_li(4,"1/2"), my_li(4,"1/9")
Li4m2 = my_li4_neg2()
Li3mh = my_li(3,"-1/2")

# ID2: R3^(2) + zeta(2)/2 = -Li2(1/3) - Li2(-1/3) + pi^2/12
id2 = R32 + z2/2 - (-Li21 - Li2m + pi**2/12)
log(f"ID2 residual (250 dps, my HPLs): {mp.nstr(abs(id2),3)}")
# ID3: R3^(3) + zeta(2) Psi^(1) - zeta(3)/2 = W3closed
W3cl = (-26*z3 - 87*Li31 - 130*Li3m + mp.mpf(21)/4*pi**2*L2 - mp.mpf(43)/12*pi**2*L3
        + 67*Li21*L2 - 45*Li21*L3 + mp.mpf(65)/3*L2**3 - 65*L2**2*L3 + mp.mpf(109)/2*L2*L3**2
        - mp.mpf(37)/3*L3**3)
id3 = R33 + z2*P1 - z3/2 - W3cl
log(f"ID3 residual: {mp.nstr(abs(id3),3)}")
# ID4: R3^(4) + zeta(2) Psi^(2) + (2/3) zeta(3) log2 = W4closed (= F4 of COR7)
W4cl = (-mp.mpf(22)/3*Li4m2 - mp.mpf(17)/2*Li41 + mp.mpf(20)/3*Li423 + mp.mpf(110)/9*Li412
        + mp.mpf(25)/48*Li419 - 12*Li31*L2 + mp.mpf(28)/3*Li31*L3 - mp.mpf(16)/3*Li3m*L2
        + mp.mpf(16)/3*Li3m*L3 + mp.mpf(14)/3*z3*L2 - mp.mpf(25)/9*z3*L3 - 4*Li21**2
        + mp.mpf(13)/3*Li21*L2**2 - 6*Li21*L2*L3 + mp.mpf(2)/3*Li21*L3**2
        + mp.mpf(11)/18*pi**2*Li21 - mp.mpf(31)/27*pi**2*L2**2 + mp.mpf(31)/72*pi**2*L3**2
        - mp.mpf(271)/1296*pi**4 + mp.mpf(46)/27*L2**4 - mp.mpf(62)/9*L2**3*L3
        + mp.mpf(61)/6*L2**2*L3**2 - mp.mpf(17)/3*L2*L3**3 + mp.mpf(103)/144*L3**4)
id4 = R34 + z2*P2 + mp.mpf(2)/3*z3*L2 - W4cl
log(f"ID4 residual: {mp.nstr(abs(id4),3)}")
assert abs(id2) < mp.mpf('1e-225') and abs(id3) < mp.mpf('1e-225') and abs(id4) < mp.mpf('1e-225')
with open('progress.log','a') as f:
    f.write("referee step(ID 250dps): ID2/3/4 residuals %.1e/%.1e/%.1e with MY hpl_value(.,1/3) + MY Li values\n" % (abs(id2),abs(id3),abs(id4)))
log("STEP PASS")
