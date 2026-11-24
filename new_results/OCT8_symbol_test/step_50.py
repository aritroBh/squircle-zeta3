import mpmath as mp
mp.mp.dps = 100
for x in [3, 80, 5, 125/3]:
    lhs = mp.polylog(3, -x) - mp.polylog(3, -1/x)
    rhs = -(mp.log(x)**3 + mp.pi**2*mp.log(x))/6
    print(x, "inversion residual:", mp.nstr(abs(lhs-rhs), 3))
