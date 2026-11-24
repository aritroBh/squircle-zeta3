"""step_26.py -- PART 3b: numeric check M2 t-integral vs corner_series M2 (3 angles, 60 digits)."""
import mpmath as mp
from corner_series import corner_coeffs
mp.mp.dps = 80

def M2_t_integral(a):
    c = mp.cos(a/2); s = mp.sin(a/2); ca = mp.cos(a)
    r1, r2, t1 = 2*c - 1, -1 - 2*c, 1 + 2*s
    def P(t): return t*t - 2*ca*t + 1
    def RA2(t):
        # use the derived factored form: build from definition for safety
        g = (t*t - 1)/(2*(t - ca))
        S0 = P(t)**2/(4*(t - ca)**2)
        J = P(t)/((t*t - 1)*(t - ca))
        A0 = P(t)/((t - r1)*(t - r2))
        cag_g2 = -(t*t - 1)*(t*t - 2*ca*t + (2*ca*ca - 1))/(4*(t - ca)**2)
        g2_2_cag = (t*t - 1)*(t*t - ca*t + ca*ca - 1)/(2*(t - ca)**2)
        R1c = cag_g2/S0
        R2c = g2_2_cag/(2*S0) - cag_g2**2/(2*S0**2)
        return J*(A0*(R2c + R1c + mp.mpf('0.5')) - mp.mpf('0.5'))
    def RB2(t):
        g = (t*t - 1)/(2*(t - ca))
        S0 = P(t)**2/(4*(t - ca)**2)
        J = P(t)/((t*t - 1)*(t - ca))
        A0 = P(t)/((t - r1)*(t - r2))
        cag_g2 = -(t*t - 1)*(t*t - 2*ca*t + (2*ca*ca - 1))/(4*(t - ca)**2)
        return -J*A0*(1 + cag_g2/S0)
    def RC2(t):
        g = (t*t - 1)/(2*(t - ca))
        J = P(t)/((t*t - 1)*(t - ca))
        A0 = P(t)/((t - r1)*(t - r2))
        return J*A0/2
    def integrand(t):
        g = (t*t - 1)/(2*(t - ca))
        lg = mp.log(g); l1g = mp.log(1 + g)
        return RA2(t)*lg**2 + RB2(t)*lg*l1g + RC2(t)*l1g**2
    pts = [1, 1 + (t1-1)/64, 1 + (t1-1)/8, 1 + (t1-1)/2, t1]
    return mp.quad(integrand, pts)

for nm, a in [("3pi/4", 3*mp.pi/4), ("pi/2", mp.pi/2), ("1.3", mp.mpf('1.3'))]:
    v = M2_t_integral(a)
    ref = corner_coeffs(a, 2)[2]
    print(f"  {nm:>6}: t-integral - series M2 = {mp.nstr(v - ref, 3)}  ({float(-mp.log10(abs(v-ref))):.1f} dig)")
