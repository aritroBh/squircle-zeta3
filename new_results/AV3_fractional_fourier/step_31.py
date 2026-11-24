from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi)
mp.dps = 80
def N(x,d=40): return nstr(x,d)
Lp = log(pi)
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*Lp*(-1)**n/n + (2*Jc + Kc)/n )
# manual partial sums with known asymptotic tail: bp2(2k+1) ~ 2 pi Lp/(2k+1) + O(log k / k^3)
M = 20000
part = nsum(lambda k: bp2(2*k+1)*(-1)**k, [0, M])
# tail: alternating a(k)=bp2(2k+1): leading 2πLp/(2k+1): tail sum ≈ (-1)^{M+1} * a(M+1)/2 approx
print("partial to M=", M, ":", N(part,30))
est_tail = (-1)**(M+1) * bp2(2*(M+1)+1)/2
print("tail est:", nstr(est_tail,8), " => total ~", N(part+est_tail,30))
# nsum on the piece-sum again but in two halves to cross-check
half1 = nsum(lambda k: bp2(2*k+1)*(-1)**k, [0, 300])
half2 = nsum(lambda k: bp2(2*k+1)*(-1)**k, [301, inf])
print("nsum total:", N(half1+half2,30))
print("exact     :", N((pi/2)**2*log(pi/2),30))
# and the P-pieces total
beta = lambda s: (zeta(s, mpf(1)/4) - zeta(s, mpf(3)/4))/4**s
def betap(s):  return -nsum(lambda k: (-1)**k*log(2*k+1)/(2*k+1)**s, [0, inf])
Cisum = nsum(lambda k: (-1)**k*ci((2*k+1)*pi)/(2*k+1)**3, [0, inf])
tot = pi**2*Lp/2 -(4/pi)*(Cisum + (Lp-euler)*beta(3) + betap(3)) -(8/pi)*beta(4) -pi**2/8
print("P-pieces  :", N(tot,30))
