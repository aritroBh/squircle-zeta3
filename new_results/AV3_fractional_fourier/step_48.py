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
Pp3 = -(pi**3/2)*(Lp+mpf(1)/3) + (12/pi)*(CiS4 + euler*eta4 - Lp*zeta(4) + etap4) - (mpf(165)/(4*pi))*zeta(4)
print("P'(3) assembled =", N(Pp3,40))
print("P'(3) exact     =", N(Pp3_exact,40), " diff:", nstr(abs(Pp3-Pp3_exact),3))
# CiS4 closed form solved from exact P'(3):
CiS4_closed = (pi/12)*( Pp3_exact + (pi**3/2)*(Lp+mpf(1)/3) + (mpf(165)/(4*pi))*zeta(4) ) - euler*eta4 + Lp*zeta(4) - etap4
print("CiS4 = sum (-1)^n Ci(n pi)/n^4")
print("direct =", N(CiS4,42))
print("closed =", N(CiS4_closed,42), " diff:", nstr(abs(CiS4-CiS4_closed),3))
# simplify: symbolic
import sympy as sp
z4, zp4, g_, L_, l2_ = sp.symbols('z4 zp4 gam L l2')
expr = (sp.pi/12)*( sp.pi**3*(12*L_+1)/16 + (sp.pi**3/2)*(L_+sp.Rational(1,3)) + sp.Rational(165,4)*z4/sp.pi ) - g_*sp.Rational(7,8)*z4 + L_*z4 - (l2_*z4/8 + sp.Rational(7,8)*zp4)
print("CiS4 =", sp.simplify(expr))
logline("step_43: CiS4 closed form via P'(3), involves zeta'(4); verified")
