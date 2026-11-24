from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi, quad)
mp.dps = 100
def N(x,d=40): return nstr(x,d)
beta = lambda s: (zeta(s, mpf(1)/4) - zeta(s, mpf(3)/4))/4**s
def betap(s):  return -nsum(lambda k: (-1)**k*log(2*k+1)/(2*k+1)**s, [0, inf])
Lp = log(pi)
b3, b4 = beta(3), beta(4); bp3 = betap(3)
Cisum = nsum(lambda k: (-1)**k*ci((2*k+1)*pi)/(2*k+1)**3, [0, inf])
P1 = pi**2*Lp/2
P2 = -(4/pi)*(Cisum + (Lp-euler)*b3 + bp3)
P3 = -(8/pi)*b4
P4 = -pi**2/8
for nm,v in [("P1",P1),("P2",P2),("P3",P3),("P4",P4),("Cisum",Cisum),("bp3",bp3),("b4",b4)]:
    print(nm, N(v,36))
print("sum =", N(P1+P2+P3+P4,36), " exact =", N((pi/2)**2*log(pi/2),36))
# piece-check P3,P4,P1 against direct sums
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*Lp*(-1)**n/n + (2*Jc + Kc)/n )
chkP4 = nsum(lambda k: (2/pi)*(-2)/(2*k+1)**3 * (-1)**k, [0, inf])
chkP3 = nsum(lambda k: (2/pi)*(-4)/(2*k+1)**4 * (-1)**k, [0, inf])
chkP1 = nsum(lambda k: (2/pi)*pi**2*Lp/(2*k+1) * (-1)**k, [0, inf])
print("chkP4", N(chkP4,30), " chkP3", N(chkP3,30), " chkP1", N(chkP1,30))
chkP2 = nsum(lambda k: (2/pi)*(-2)*(ci((2*k+1)*pi)-euler+Lp-log(2*k+1))/(2*k+1)**3 * (-1)**k, [0, inf])
print("chkP2", N(chkP2,36), " P2 =", N(P2,36), " diff:", nstr(abs(chkP2-P2),3))
