# step_52: golden chi3 identity + octagon D8/L8, two precisions
import mpmath as mp
def setup(dps):
    mp.mp.dps = dps
    def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
    def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
    def L3(x):
        x = mp.mpf(x)
        return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
    L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
    L8 = (mp.zeta(3, mp.mpf(1)/8) - mp.zeta(3, mp.mpf(3)/8) - mp.zeta(3, mp.mpf(5)/8) + mp.zeta(3, mp.mpf(7)/8))/512
    return Li3, Li2, L3, L5, L8
for dps in (150, 250):
    Li3, Li2, L3, L5, L8 = setup(dps)
    phi = (1+mp.sqrt(5))/2; lphi = mp.log(phi)
    chi3g = (Li3(1/phi)-Li3(-1/phi))/2
    # identity: chi3(1/phi) = (25 sqrt5 / 48) L5 - pi^2 lphi/12 + lphi^3/4
    rhs = mp.sqrt(5)*mp.mpf(25)/48*L5 - mp.pi**2*lphi/12 + lphi**3/4
    print(f"dps={dps}: golden chi3 identity residual:", abs(chi3g - rhs))
    b = mp.sqrt(2)-1
    D8 = L3(b)-L3(-b)
    print(f"dps={dps}: D8/L8 =", mp.nstr(D8/L8, 40))
    import sympy as sp
    print("   nsimplify:", sp.nsimplify(str(D8/L8), [sp.sqrt(2)]))
    # golden D/L5 = 25 sqrt5/24 check
    D = L3(1/phi)-L3(-1/phi)
    print("   D/L5 - 25sqrt5/24:", D/L5 - 25*mp.sqrt(5)/24)
