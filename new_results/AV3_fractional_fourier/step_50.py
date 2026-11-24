from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler)
mp.dps = 100
Lp = log(pi)
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*Lp*(-1)**n/n + (2*Jc + Kc)/n )
def b_c2(n): return 2*pi*(-1)**(n+1)/n - 4*(1-(-1)**n)/(pi*n**3)
def ap3(n): return -(3*bp2(n) + b_c2(n))/n
def ap3_pred(n):
    return ( 6*pi*Lp*(-1)**n/n**2 + 2*pi*(-1)**n/n**2
             + (12/pi)*(ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n**4
             + (22/pi)*(1-(-1)**n)/n**4 )
for n in [1,2,3,4]:
    print(n, nstr(ap3(n),22), nstr(ap3_pred(n),22), nstr(abs(ap3(n)-ap3_pred(n)),3))
