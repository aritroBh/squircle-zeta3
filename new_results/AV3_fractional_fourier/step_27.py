from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi)
mp.dps = 100
def N(x,d=40): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
beta = lambda s: (zeta(s, mpf(1)/4) - zeta(s, mpf(3)/4))/4**s
def betap(s):  return -nsum(lambda k: (-1)**k*log(2*k+1)/(2*k+1)**s, [0, inf])
G = beta(2); Lp = log(pi)

gfun = lambda n: (-1)**(n+1)*(si(n*pi) - pi/2)   # >0, smooth
resid = nsum(lambda k: gfun(2*k+1)*(-1)**k/(2*k+1)**2, [0, inf])
S_Sisin = (pi/2)*G + resid
print("S_Sisin =", N(S_Sisin,42))
print("(pi^2/4)log2 =", N(pi**2/4*log(2),42), " diff:", nstr(abs(S_Sisin-pi**2/4*log(2)),3))

# T'(2):  T'(2) = (pi^2/2)log pi - (4/pi)[Cisum + (log pi - gamma) beta(3) + beta'(3)] - (8/pi)beta(4) - pi^2/8
Cisum = nsum(lambda k: (-1)**k*ci((2*k+1)*pi)/(2*k+1)**3, [0, inf])
b3 = beta(3); b4 = beta(4); bp3 = betap(3)
Tp2_series = (pi**2/2)*Lp - (4/pi)*(Cisum + (Lp-euler)*b3 + bp3) - (8/pi)*b4 - pi**2/8
Tp2_exact  = (pi/2)**2*log(pi/2)
print("T'(2) assembled =", N(Tp2_series,40)); print("T'(2) exact     =", N(Tp2_exact,40), " diff:", nstr(abs(Tp2_series-Tp2_exact),3))
# hence Cisum closed form (solved from exact T'(2)):
Cisum_closed = (pi/4)*( (pi**2/2)*Lp - pi**2/8 - (8/pi)*b4 - Tp2_exact ) - (Lp-euler)*b3 - bp3
print("Cisum direct =", N(Cisum,42))
print("Cisum closed =", N(Cisum_closed,42), " diff:", nstr(abs(Cisum-Cisum_closed),3))
logline("step_26: S_Sisin=(pi^2/4)log2 and T'(2) assembly with beta'(3), beta(4), Cisum closed form verified")
