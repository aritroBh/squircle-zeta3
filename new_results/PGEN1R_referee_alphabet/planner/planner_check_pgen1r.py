# Planner check of PGEN1R's identities C and D (own code; L(chi_5,3) from Hurwitz zeta).
import mpmath as mp
mp.mp.dps = 90
s5 = mp.sqrt(5); phi = (1 + s5)/2; r = 1/phi; lp = mp.log(phi)
chi = {1: 1, 2: -1, 3: -1, 4: 1}
L5 = sum(chi[a]*mp.zeta(3, mp.mpf(a)/5) for a in chi)/mp.mpf(5)**3
chi3 = lambda x: (mp.polylog(3, x) - mp.polylog(3, -x))/2
D = chi3(r) - ((25*s5/48)*L5 - (mp.pi**2/12)*lp + lp**3/4)
print("D: chi3(1/phi) - [(25 sqrt5/48) L(chi5,3) - (pi^2/12) log phi + (1/4) log^3 phi] =", mp.nstr(D, 5))
def L3(x):  # Zagier's single-valued trilogarithm, real x
    lx = mp.log(abs(x))
    return mp.re(mp.polylog(3, x) - lx*mp.polylog(2, x) + lx**2*(-mp.log(1 - x))/3)
c1 = -mp.mpf(75)/32 + 81*s5/64; c2 = -mp.mpf(35)/32 + 53*s5/64
c3 = mp.mpf(33)/64 - 19*s5/80; c3loc = mp.mpf(5)/32 - s5/8
N3 = c1*L3(r) + c2*L3(-r) + (c3 - c3loc)*mp.zeta(3)
a = mp.mpf(1)/64 + 31*s5/320; b = mp.mpf(875)/768 - 125*s5/192
print("C: N3 - a zeta(3) - b L(chi5,3) =", mp.nstr(N3 - a*mp.zeta(3) - b*L5, 5), "  (N3 =", mp.nstr(N3, 20), ")")
