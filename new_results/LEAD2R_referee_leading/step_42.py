import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 60
from lead2r_lib import Efun
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
s = mp.mpf('0.99')
# A: direct on [0,1] with strong subtraction of 2 Taylor terms: h(x) = 1 - x/2 - x^2/2 + ...
def I2(t, f, c1, c2):
    # int_0^1 x^t f(x) dx, f = 1 + c1 x + c2 x^2 + O(x^3)
    return 1/(t+1) + c1/(t+2) + c2/(t+3) + mp.quad(lambda x: x**t*(f(x)-1-c1*x-c2*x**2), [0, mp.mpf('0.5'), 1])
# h(x) = ((1+x)/sqrt(1+x^2)-1)/x: expand: (1+x)(1 - x^2/2 + 3x^4/8) = 1 + x - x^2/2 - x^3/2 + 3x^4/8; -1; /x = 1 - x/2 - x^2/2 + 3x^3/8
F1 = I2(-s, h, mp.mpf('-0.5'), mp.mpf('-0.5'))
F2 = I2(s, h, mp.mpf('-0.5'), mp.mpf('-0.5'))
print("F1+F2 (Taylor-subtracted):", mp.nstr(F1+F2, 15))
rhs = (Efun(s)-Efun(-s))/s
print("rhs:", mp.nstr(rhs, 15))
print("err:", mp.nstr(abs(F1+F2-rhs),3))
# B: e^-r approach, print both pieces
J1 = 1/(1-s) + mp.quad(lambda r: mp.exp(-(1-s)*r)*(h(mp.exp(-r))-1), [0, 20, mp.inf])
J2 = 1/(1+s) + mp.quad(lambda r: mp.exp(-(1+s)*r)*(h(mp.exp(-r))-1), [0, 20, mp.inf])
print("J1+J2:", mp.nstr(J1+J2,15), "err:", mp.nstr(abs(J1+J2-rhs),3))
