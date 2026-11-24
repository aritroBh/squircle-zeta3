import sys
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from qsqrt2 import QS2 as Q
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 50

# s = sin(pi/8), c = cos(pi/8) for m=4 (regular octagon, 2m-gon m=4)
# s^2 = (2-sqrt2)/4 exactly; c = u/2 with u = sqrt(2+sqrt2)
R2 = Q(0,1)
s2 = (2-R2)/4
s4 = s2*s2
s6 = s4*s2
# C_m^(5) = 2880 m s^4 c (65 s^6 - 63 s^4 + 17 s^2 - 1), m=4
A5 = 65*s6 - 63*s4 + 17*s2 - 1
# C_m^(2,3) = -960 m c s^4 (45 s^6 - 27 s^4 + 6 s^2 - 1)
A23 = 45*s6 - 27*s4 + 6*s2 - 1
print("65s^6-63s^4+17s^2-1 =", A5)
print("45s^6-27s^4+6s^2-1 =", A23)
# C5 = 2880*4*s4*A5 * c  = (coeff in Q(sqrt2)) * c ; c = u/2, u = sqrt(2+sqrt2)
C5_body = 2880*4*s4*A5      # times c
C23_body = -960*4*s4*A23    # times c
print("C_4^(5)   = (%s) * c" % C5_body)
print("C_4^(2,3) = (%s) * c" % C23_body)
u = mp.sqrt(2+mp.sqrt(2)); c = u/2
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q))*mp.sqrt(2)
print("numeric C_4^(5)   =", mp.nstr(q2f(C5_body)*c, 30))
print("numeric C_4^(2,3) =", mp.nstr(q2f(C23_body)*c, 30))
# express via u: c = u/2 -> C5 = (C5_body/2)*u
print("C_4^(5)   = (%s) * u" % (C5_body/2))
print("C_4^(2,3) = (%s) * u" % (C23_body/2))
