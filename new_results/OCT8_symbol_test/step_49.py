import mpmath as mp
mp.mp.dps = 100
# inversion for Li_3 at negative argument: Li_3(-x) + Li_3(-1/x) = -(1/6)ln^3(x) - (pi^2/6)ln x + C?
# determine C at x=1: 2 Li_3(-1) = 2*(-3 zeta3/4) = -3/2 zeta3 ; RHS at x=1: C => C = -3/2 zeta(3)?? 
for x in [3, 80, 5]:
    lhs = mp.polylog(3, -x) + mp.polylog(3, -1/x)
    rhs = -(mp.log(x)**3)/6 - (mp.pi**2)*mp.log(x)/6 - 3*mp.zeta(3)/2
    print(x, mp.nstr(abs(lhs-rhs), 3))
