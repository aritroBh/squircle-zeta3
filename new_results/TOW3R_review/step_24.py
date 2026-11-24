import sympy as sp
d, m = sp.symbols('d m'); i = sp.I
N5 = ( -16*i*d**7 + 304*i*d**6 + (-80-2520*i)*d**5 + 160*d**4*m + (1040+11840*i)*d**4
      - 2080*d**3*m + (-5140-34069*i)*d**3 + (10280-240*i)*d**2*m + (11440+60081*i)*d**2
      + 240*i*d*m**2 + (-22880+960*i)*d*m + (-9600-59920*i)*d - 960*i*m**2 + 19200*m + 25200*i )
a5 = N5/(480*(2*d-5)*(2*d-3))
v = sp.expand(a5.subs({d:5,m:3}))
print('a5(5,3) =', v)
from fractions import Fraction as Fr
ratio = (Fr(3,56), Fr(-31,70))
print('ratio data =', ratio, ' 2*a5 =', 2*v)
print('equal to a5?', sp.Rational(ratio[0])+sp.I*ratio[1] == v)
print('equal to 2*a5?', sp.Rational(ratio[0])+sp.I*ratio[1] == 2*v)
