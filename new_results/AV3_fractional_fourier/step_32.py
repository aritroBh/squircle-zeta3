from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi)
mp.dps = 80
Lp = log(pi)
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*Lp*(-1)**n/n + (2*Jc + Kc)/n )
def raw_odd(n):
    return (2/pi)*( pi**2*Lp/n - 2*(ci(n*pi)-euler+Lp-log(n))/n**3 - 4/n**4 - 2/n**3 )
for k in range(6):
    n = 2*k+1
    print(n, nstr(bp2(n),20), nstr(raw_odd(n),20), nstr(abs(bp2(n)-raw_odd(n)),3))
