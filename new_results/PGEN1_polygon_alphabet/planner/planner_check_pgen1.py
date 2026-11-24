# Planner check of PGEN1: (1) the golden chi_3 identity, (2) the eq:chi3 control, (3) L3(1/2) = 7/8 zeta(3) with Zagier's L3.
import mpmath as mp
mp.mp.dps = 120
s5 = mp.sqrt(5); phi = (1 + s5)/2; L = mp.log(phi)
def Lchi(D, s, chi):
    return mp.mpf(D)**(-s) * mp.fsum(chi(a) * mp.zeta(s, mp.mpf(a)/D) for a in range(1, D))
chi5 = lambda a: [0, 1, -1, -1, 1][a % 5]
chi8 = lambda a: {1: 1, 7: 1, 3: -1, 5: -1}.get(a % 8, 0)
L5 = Lchi(5, 3, chi5); L8 = Lchi(8, 3, chi8)
direct5 = mp.nsum(lambda n: chi5(int(n)) / n**3, [1, mp.inf]) if False else None
chi3 = lambda x: (mp.polylog(3, x) - mp.polylog(3, -x)) / 2
lhs = chi3(1/phi); rhs = 25*s5/48*L5 - mp.pi**2/12*L + L**3/4
print("L(chi5,3) =", mp.nstr(L5, 25), "  partial-sum check:", mp.nstr(abs(L5 - mp.fsum(chi5(n)/mp.mpf(n)**3 for n in range(1, 200001))), 3))
print("golden chi3 identity: |chi3(1/phi) - [(25 sqrt5/48) L(chi5,3) - pi^2/12 log phi + log^3 phi/4]| =", mp.nstr(abs(lhs - rhs), 3))
print("L(chi8,3) =", mp.nstr(L8, 25), " (PGEN1 prints 0.958380454563094562...)")
L3 = lambda x: mp.re(mp.polylog(3, x)) - mp.log(abs(x))*mp.re(mp.polylog(2, x)) - mp.log(abs(x))**2*mp.log(abs(1 - x))/3
print("Zagier L3(1/2) - 7/8 zeta(3) =", mp.nstr(L3(mp.mpf(1)/2) - mp.mpf(7)/8*mp.zeta(3), 3))
