
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 1200
# quick reuse of step_08 machinery
exec(open("steps/step_08.py").read().split("DPS = 1200")[0])
def my_polylog_series(s, z):
    z = mp.mpc(z); az = abs(z)
    eps = mp.mpf(10) ** (-1230)
    tot = mp.mpc(0); n = 1; zn = z
    while True:
        t = zn / mp.mpf(n) ** s
        tot += t
        if n > 100 and abs(t) * az / (1 - az) < eps: break
        n += 1; zn *= z
    return tot
def bernoulli4(y): return y**4 - 2*y**3 + y**2 - mp.mpf(1)/30
def my_li4_neg2():
    y = mp.mpc(mp.mpf(1)/2) + mp.log(2) / (2 * mp.pi * mp.j)
    return -(2 * mp.pi * mp.j) ** 4 / mp.factorial(4) * bernoulli4(y) - my_polylog_series(4, mp.mpf(-1)/2)
pi, z2, z3 = mp.pi, mp.zeta(2), mp.zeta(3)
L2, L3 = mp.log(2), mp.log(3)
Li21 = my_polylog_series(2, mp.mpf(1)/3)
Li31 = my_polylog_series(3, mp.mpf(1)/3); Li3m = my_polylog_series(3, mp.mpf(-1)/2)
Li41 = my_polylog_series(4, mp.mpf(1)/3); Li4m3 = my_polylog_series(4, mp.mpf(-1)/3)
Li423 = my_polylog_series(4, mp.mpf(2)/3); Li412 = my_polylog_series(4, mp.mpf(1)/2)
Li419 = my_polylog_series(4, mp.mpf(1)/9); Li4m2 = my_li4_neg2()
W2 = L2**2/2 + Li21 + mp.log(mp.mpf(3)/2)**2/2
W3 = mp.mpf(-15)/4*z3 - 3*Li31 - 8*Li3m - pi**2*L2/12 + 4*Li21*L2 - mp.mpf(3)/2*Li21*L3 \
     + L2**3/3 - mp.mpf(5)/2*L2**2*L3 + 2*L2*L3**2 - L3**3/4
F4 = (-21*Li4m2 - mp.mpf(65)/8*Li41 + mp.mpf(13)/2*Li423 + mp.mpf(137)/6*Li412 + mp.mpf(39)/64*Li419
      + mp.mpf(15)/2*Li31*L2 + mp.mpf(37)/2*Li3m*L2 - mp.mpf(5)/4*Li3m*L3
      + mp.mpf(119)/16*z3*L2 - mp.mpf(37)/32*z3*L3
      + Li21**2 - mp.mpf(23)/4*Li21*L2**2 + mp.mpf(7)/4*Li21*L2*L3 + Li21*L3**2 - pi**2*Li21/6
      - mp.mpf(455)/144*pi**2*L2**2 + mp.mpf(9)/8*pi**2*L2*L3 - mp.mpf(7)/32*pi**2*L3**2
      - mp.mpf(949)/2160*pi**4 - mp.mpf(211)/144*L2**4 + mp.mpf(101)/24*L2**3*L3
      - mp.mpf(15)/8*L2**2*L3**2 - mp.mpf(35)/24*L2*L3**3 + mp.mpf(113)/192*L3**4)
psi2_closed = mp.mpf(3)/32*pi**2 - mp.mpf(9)/8*Li21 + mp.mpf(3)/4*L2*L3 - mp.mpf(9)/16*L3**2
W4 = F4 + z2 * psi2_closed + mp.mpf(3)/8 * z3 * L3
M3 = W4 + W3 + W2
m3ref = open("../COR5_corner_classicality/M3_2pi3_1100.txt").read().strip()
diff = abs(M3 - mp.mpf(m3ref))
print("fully-explicit M3 (certified Psi2) vs COR5 1095-digit value:", mp.nstr(diff, 3))
print("digits of agreement:", int(-mp.log10(diff)))
