
from fractions import Fraction as F
import dec5d_lib as D
QS5 = D.QS5; PHI = D.PHI
s = 1/(2*PHI)                     # sin(pi/10) in Q(sqrt5)
cOverL = PHI*(3*D.SQ5+5)/4        # cos(pi/10)/L in Q(sqrt5)
poly5 = 65*s**6 - 63*s**4 + 17*s**2 - 1
poly23 = 45*s**6 - 27*s**4 + 6*s**2 - 1
C5_5_over_L = QS5(2880*5)*s**4*poly5*cOverL
C5_23_over_L = QS5(-960*5)*s**4*poly23*cOverL
print("local-block zeta(5) coeff in P5/L   =", C5_5_over_L)
print("local-block zeta(2)zeta(3) coeff in P5/L =", C5_23_over_L)
import mpmath as mp
mp.mp.dps = 40
v1 = float(C5_5_over_L); v2 = float(C5_23_over_L)
print("numeric:", mp.nstr(v1, 30), "vs 65.59694677451535325292533845872515691473")
print("numeric:", mp.nstr(v2, 30), "vs 131.4394994149313640595416717604309158838")
