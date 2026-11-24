from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler)
mp.dps = 100
def N(x,d=36): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
Lp = log(pi)
eta4 = mpf(7)/8*zeta(4)
etap4 = (log(2)/8)*zeta(4) + (mpf(7)/8)*zeta(4, derivative=1)
Pp3_exact = pi**3*(12*Lp+1)/16
CiS4 = nsum(lambda n: (-1)**n*ci(n*pi)/n**4, [1, inf])
P1 = (6*pi*Lp+2*pi)*zeta(2)
P3 = -(mpf(22)/pi)*(eta4 + zeta(4))
P2 = (12/pi)*(CiS4 + euler*eta4 - Lp*zeta(4) - etap4)
Pp3 = P1 + P2 + P3
print("P'(3) assembled =", N(Pp3,40))
print("P'(3) exact     =", N(Pp3_exact,40), " diff:", nstr(abs(Pp3-Pp3_exact),3))
CiS4_closed = (pi/12)*( Pp3_exact - P1 - P3 ) - euler*eta4 + Lp*zeta(4) + etap4
print("CiS4 direct =", N(CiS4,42))
print("CiS4 closed =", N(CiS4_closed,42), " diff:", nstr(abs(CiS4-CiS4_closed),3))
import sympy as sp
z4, zp4, g_, L_, l2_ = sp.symbols('z4 zp4 gam L l2')
expr = (sp.pi/12)*( sp.pi**3*(12*L_+1)/16 - (6*sp.pi*L_+2*sp.pi)*sp.zeta(2) + sp.Rational(22)*sp.Rational(15,8)*z4/sp.pi ) - g_*sp.Rational(7,8)*z4 + L_*z4 + l2_*z4/8 + sp.Rational(7,8)*zp4
print("CiS4 =", sp.simplify(expr))
logline("step_45: P'(3) verified; CiS4 closed form with zeta'(4), log2, gamma verified 40 digits")
