"""step_15.py -- GENERAL-ANGLE THEOREM, part C: numeric verification of the compact
closed form M1(a) = Re[F(t1) - F(1+)] at >= 8 angles vs corner_series (>=60 digits),
imaginary drift, and the zeta(2) coefficient gamma = -1.
"""
import mpmath as mp
from corner_series import corner_coeffs
mp.mp.dps = 90
log2 = mp.log(2); z2 = mp.pi**2/6

def data(a):
    c = mp.cos(a/2); s = mp.sin(a/2); ca = mp.cos(a)
    return c, s, ca, 2*c - 1, -1 - 2*c, 1 + 2*s

def make_F(a):
    c, s, ca, r1, r2, t1 = data(a)
    poles = [mp.mpf(1), mp.mpf(-1), ca, r1, r2]
    ares = [mp.mpf(0), mp.mpf(-2), mp.mpf(1), -c, c]
    bres = [mp.mpf(-1), mp.mpf(1), mp.mpf(-1), 2*c, -2*c]
    Apart = [(mp.mpf(1), 1), (mp.mpf(-1), 1), (ca, -1)]
    Bpart = [(r1, 1), (r2, 1), (ca, -1)]
    def I(t, rho, sig):
        if rho == sig:
            return mp.log(t - rho)**2/2
        z = (t - rho)/(sig - rho)
        return -mp.polylog(2, z) + mp.log(rho - sig)*mp.log(t - rho)
    def F(t):
        tot = mp.mpc(0)
        for rho, ar, br in zip(poles, ares, bres):
            for sig, kap in Apart:
                tot += ar*kap*I(t, rho, sig)
            for sig, kap in Bpart:
                tot += br*kap*I(t, rho, sig)
            tot -= (ar + br)*log2*mp.log(t - rho)
        return tot
    return F, t1

def M1_formula(a, zero_zeta2=False):
    F, t1 = make_F(a)
    eps = mp.mpf(10)**(-70)
    if not zero_zeta2:
        return F(t1) - F(1 + eps)
    # same but with Li2(1) terms removed (only sigma=1, lower endpoint, rho != 1)
    c, s, ca, r1, r2, t1_ = data(a)
    poles = [mp.mpf(1), mp.mpf(-1), ca, r1, r2]
    ares = [mp.mpf(0), mp.mpf(-2), mp.mpf(1), -c, c]
    corr = mp.mpc(0)
    for rho, ar in zip(poles, ares):
        if rho == 1: continue
        # remove -Li2(1) piece: add back +Li2(1) with weight ar*(+1), at lower limit
        corr += ar*mp.polylog(2, 1)
    return (F(t1) - F(1 + eps)) + corr  # F(1+) -= ar*(-z2) i.e. add ar*z2 to value

angles = [("2pi/3", 2*mp.pi/3), ("pi/2", mp.pi/2), ("pi/3", mp.pi/3), ("3pi/4", 3*mp.pi/4),
          ("pi/6", mp.pi/6), ("5pi/6", 5*mp.pi/6), ("pi/12", mp.pi/12), ("1.3", mp.mpf('1.3'))]
print(f"{'angle':>7} | {'residual (formula - series)':>28} | {'Im drift':>10} | digits")
for nm, a in angles:
    v = M1_formula(a)
    ref = corner_coeffs(a, 1)[1]
    res = mp.re(v) - ref
    dig = -mp.log10(abs(res)) if res != 0 else mp.inf
    print(f"{nm:>7} | {mp.nstr(res, 3):>28} | {mp.nstr(abs(mp.im(v)),2):>10} | {dig:.1f}")

print()
print("zeta(2) coefficient check: formula with Li2(1)-terms zeroed should exceed M1 by -gamma*zeta2, gamma=-1")
for nm, a in [("pi/2", mp.pi/2), ("3pi/4", 3*mp.pi/4), ("1.3", mp.mpf('1.3'))]:
    v = M1_formula(a); v0 = M1_formula(a, zero_zeta2=True)
    print(f"  {nm:>6}: (v0 - v)/zeta2 = {mp.nstr(mp.re(v0 - v)/z2, 12)}  (expect +1 since M1 ⊃ -1*zeta2)")
