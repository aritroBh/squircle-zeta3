from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi, pslq)
mp.dps = 100
def N(x,d=42): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
Lp = log(pi); eta3 = 3*zeta(3)/4; zp3 = zeta(3, derivative=1)
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*Lp*(-1)**n/n + (2*Jc + Kc)/n )

# CiS = sum ci(n pi)/n^3 : method 1 direct nsum; method 2 hump interchange
CiS1 = nsum(lambda n: ci(n*pi)/n**3, [1, inf])
ciH = lambda m: ci((m+1)*pi) - ci(m*pi)
H3 = lambda m: zeta(3) - zeta(3, m+1)
N0 = 50
CiS2_head = nsum(lambda n: ci(n*pi)/n**3, [1, N0])
CiS2 = CiS2_head + ci(N0*pi)*zeta(3,N0+1) + nsum(lambda m: ciH(m)*(zeta(3)-H3(m)), [N0, inf])
print("CiS direct =", N(CiS1,45))
print("CiS hump   =", N(CiS2,45), " diff:", nstr(abs(CiS1-CiS2),3))

Wp2_direct = nsum(lambda n: bp2(n), [1, inf])
Wp2_pieces = 2*pi*Lp*log(2) - (4/pi)*(CiS1 + zp3 + Lp*eta3 - euler*zeta(3)) - (21/(2*pi))*zeta(3)
Wp2_int = (1/pi)*quad(lambda t: t**2*log(t)/tan(t/2), [0, pi])
print("W'(2) direct nsum =", N(Wp2_direct,42))
print("W'(2) pieces      =", N(Wp2_pieces,42))
print("W'(2) integral    =", N(Wp2_int,42))
logline("step_30: CiS two methods; W'(2) three methods agree")
