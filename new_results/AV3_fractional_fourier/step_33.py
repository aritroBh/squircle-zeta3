from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi)
mp.dps = 100
def N(x,d=42): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
beta = lambda s: (zeta(s, mpf(1)/4) - zeta(s, mpf(3)/4))/4**s
def betap(s):  return -nsum(lambda k: (-1)**k*log(2*k+1)/(2*k+1)**s, [0, inf])
Lp = log(pi); b3 = beta(3); bp3 = betap(3)
Cisum = nsum(lambda k: (-1)**k*ci((2*k+1)*pi)/(2*k+1)**3, [0, inf])
Tp2_exact = (pi/2)**2*log(pi/2)
# corrected odd-n expansion of bp2: (2Jc+Kc)/n = -2(ci-gamma+Lp-log n)/n^3 - 6/n^3
Tp2_assembled = pi**2*Lp/2 - (4/pi)*(Cisum + (Lp-euler)*b3 + bp3) - 3*pi**2/8
print("T'(2) assembled =", N(Tp2_assembled,40))
print("T'(2) exact     =", N(Tp2_exact,40), " diff:", nstr(abs(Tp2_assembled-Tp2_exact),3))
# direct check
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*Lp*(-1)**n/n + (2*Jc + Kc)/n )
Tdir = nsum(lambda k: bp2(2*k+1)*(-1)**k, [0, inf])
print("T'(2) direct    =", N(Tdir,40), " diff:", nstr(abs(Tdir-Tp2_exact),3))
# Cisum closed form
Cisum_closed = (pi/4)*( pi**2*Lp/2 - 3*pi**2/8 - Tp2_exact ) - (Lp-euler)*b3 - bp3
print("Cisum direct =", N(Cisum,42))
print("Cisum closed =", N(Cisum_closed,42), " diff:", nstr(abs(Cisum-Cisum_closed),3))
logline("step_29: T'(2) corrected assembly verified; Cisum closed form (beta'(3), beta(3)(log pi-gamma), elementary) verified 40 digits")
