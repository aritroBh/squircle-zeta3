import mpmath as mp
mp.mp.dps=60
pi=mp.pi
Lchim8=(mp.zeta(3,mp.mpf(1)/8)+mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)-mp.zeta(3,mp.mpf(7)/8))/512
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
print("Lchim8 =", mp.nstr(Lchim8,25))
print("Lchim8/pi^3 =", mp.nstr(Lchim8/pi**3,15))
print("Lchi8 =", mp.nstr(Lchi8,25))
print("Lchi8/zeta3 =", mp.nstr(Lchi8/mp.zeta(3),10))
# rational recognition
from sympy import nsimplify, Rational
print("Lchim8/pi^3 ~", nsimplify(float(Lchim8/pi**3)))
