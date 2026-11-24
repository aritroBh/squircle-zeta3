"""step_17.py -- convention comparison + corrected zeta(2) check.
D10b convention: C(rho,sig) = log(rho-sig) if rho>sig else log(sig-rho) - i*pi.
Principal convention: C = mp.log(rho-sig). Compare Im drift; verify gamma=-1 cleanly.
"""
import mpmath as mp
from corner_series import corner_coeffs
mp.mp.dps = 90
log2 = mp.log(2); z2 = mp.pi**2/6; ipi = mp.pi*1j

def build(a, conv):
    c = mp.cos(a/2); s = mp.sin(a/2); ca = mp.cos(a)
    r1, r2, t1 = 2*c - 1, -1 - 2*c, 1 + 2*s
    poles = [mp.mpf(1), mp.mpf(-1), ca, r1, r2]
    ares = [mp.mpf(0), mp.mpf(-2), mp.mpf(1), -c, c]
    bres = [mp.mpf(-1), mp.mpf(1), mp.mpf(-1), 2*c, -2*c]
    def C(rho, sig):
        if conv == "d10b":
            return mp.log(rho - sig) if rho > sig else mp.log(sig - rho) - ipi
        return mp.log(rho - sig)
    def I(t, rho, sig):
        if rho == sig:
            return mp.log(t - rho)**2/2
        return -mp.polylog(2, (t - rho)/(sig - rho)) + C(rho, sig)*mp.log(t - rho)
    def F(t):
        tot = mp.mpc(0)
        for rho, ar, br in zip(poles, ares, bres):
            tot += ar*(I(t,rho,mp.mpf(1)) + I(t,rho,mp.mpf(-1)) - I(t,rho,ca) - log2*mp.log(t-rho))
            tot += br*(I(t,rho,r1) + I(t,rho,r2) - I(t,rho,ca) - log2*mp.log(t-rho))
        return tot
    return F, t1, ares, poles

print("convention comparison (Im drift and real residual vs corner_series):")
for nm, a in [("pi/2", mp.pi/2), ("3pi/4", 3*mp.pi/4), ("pi/6", mp.pi/6), ("pi/12", mp.pi/12), ("1.3", mp.mpf('1.3'))]:
    ref = corner_coeffs(a, 1)[1]
    eps = mp.mpf(10)**(-70)
    row = f"  {nm:>6}: "
    for conv in ["d10b", "principal"]:
        F, t1, ares, poles = build(a, conv)
        v = F(t1) - F(1 + eps)
        row += f"[{conv}: resid {mp.nstr(mp.re(v)-ref,2)}, Im {mp.nstr(abs(mp.im(v)),2)}] "
    print(row)

print()
print("corrected zeta(2) check: remove Li2(1) terms => v0 = v - gamma*zeta2 with gamma=-1, i.e. (v-v0)/zeta2 = -1... ")
print("direct: weight of Li2(1) inside v = -sum_{rho!=1} a_rho = +1 in F(1+); in M1=F(t1)-F(1+): gamma = sum a_rho(rho!=1) = -1")
for nm, a in [("pi/2", mp.pi/2), ("1.3", mp.mpf('1.3'))]:
    F, t1, ares, poles = build(a, "d10b")
    eps = mp.mpf(10)**(-70)
    v = F(t1) - F(1 + eps)
    gamma = sum(ar for rho, ar in zip(poles, ares) if rho != 1)
    # zero the Li2(1) terms: in F(1+) the term is ar*(-Li2(1)); in v it is +ar*Li2(1). Remove: v0 = v - gamma*z2
    v0 = v - gamma*z2
    print(f"  {nm:>6}: (v - v0)/zeta2 = {mp.nstr(mp.re(v - v0)/z2, 10)} = gamma = {mp.nstr(gamma, 3)} (expect -1)")
