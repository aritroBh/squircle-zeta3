# 26_check_2d_theta.py -- checks used in the proof of Zb:
#  (i) 2-D theta-representation F(eps) = (1-eps)^-1 int_0^1 D2(s)(s(1-s))^(theta-1)(s^(1+theta)+(1-s)^(1+theta))^(-1-eps) ds
#  (ii) X2 = int_0^1 D2(s,1-s) log(s) ds/s = M0 - M1 ;  (iii) Zb = F(0) - X1 - X2.
import mpmath as mp, json
from squircle_consts import consts
mp.mp.dps = 40
M0, M1, M2 = consts()
D2 = lambda s: 1 - mp.sqrt(s*s + (1-s)**2)
for eps in (mp.mpf('0.1'), mp.mpf('0.37')):
    th = eps/(1-eps)
    Ft = mp.quad(lambda s: D2(s)*(s*(1-s))**(th-1)*(s**(1+th)+(1-s)**(1+th))**(-1-eps), [0, 0.5, 1])/(1-eps)
    Fd = 2*mp.quad(lambda t: t**(eps-1)+(1-t)**(eps-1)-mp.sqrt(t**(2*eps-2)+(1-t)**(2*eps-2)), [0, 0.25, 0.5])
    print("eps=%s  theta-rep %s  direct %s  diff %s" % (eps, mp.nstr(Ft, 25), mp.nstr(Fd, 25), mp.nstr(Ft-Fd, 3)))
X2 = mp.quad(lambda s: D2(s)*mp.log(s)/s, [0, 0.5, 1])
print("X2 - (M0 - M1) =", mp.nstr(X2 - (M0 - M1), 3))
X1 = mp.quad(lambda s: D2(s)*mp.log(1-s)/s, [0, 0.5, 1])
Zb = mp.mpf(json.load(open("Q1_sub_400.json"))["Zb"])
print("Zb - (F(0) - X1 - X2) =", mp.nstr(Zb - (-2*M0 - X1 - X2), 3))
