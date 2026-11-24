"""step_16.py -- rerun of step_15 with fixed formatting."""
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
    v = F(t1) - F(1 + eps)
    if zero_zeta2:
        c, s, ca, r1, r2, t1_ = data(a)
        poles = [mp.mpf(1), mp.mpf(-1), ca, r1, r2]
        ares = [mp.mpf(0), mp.mpf(-2), mp.mpf(1), -c, c]
        corr = mp.mpc(0)
        for rho, ar in zip(poles, ares):
            if rho == 1: continue
            corr += ar*mp.polylog(2, 1)
        v = v + corr
    return v

angles = [("2pi/3", 2*mp.pi/3), ("pi/2", mp.pi/2), ("pi/3", mp.pi/3), ("3pi/4", 3*mp.pi/4),
          ("pi/6", mp.pi/6), ("5pi/6", 5*mp.pi/6), ("pi/12", mp.pi/12), ("1.3", mp.mpf('1.3'))]
print(f"{'angle':>7} | {'residual':>12} | {'Im drift':>12} | digits")
for nm, a in angles:
    v = M1_formula(a)
    ref = corner_coeffs(a, 1)[1]
    res = mp.re(v) - ref
    dig = float(-mp.log10(abs(res))) if res != 0 else 999.0
    print(f"{nm:>7} | {mp.nstr(res, 3):>12} | {mp.nstr(abs(mp.im(v)),2):>12} | {dig:.1f}")

print()
print("zeta(2) coefficient: (formula with Li2(1) zeroed - formula)/zeta2 should be +1")
for nm, a in [("pi/2", mp.pi/2), ("3pi/4", 3*mp.pi/4), ("1.3", mp.mpf('1.3'))]:
    v = M1_formula(a); v0 = M1_formula(a, zero_zeta2=True)
    print(f"  {nm:>6}: {mp.nstr(mp.re(v0 - v)/z2, 12)}")
