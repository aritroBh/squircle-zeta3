# Planner check of COR5's F4 closed form at 400 digits, and of the claim that the naive-symbol Li4 combination is wrong.
import mpmath as mp
mp.mp.dps = 420
F4 = mp.mpf(open('<project>/research/08_theory/COR5_corner_classicality/hex_F4.txt').read().split()[0][:430])
L4 = lambda x: mp.polylog(4, x); L3 = lambda x: mp.polylog(3, x); L2 = lambda x: mp.polylog(2, x)
l2, l3, z3, p2 = mp.log(2), mp.log(3), mp.zeta(3), mp.pi**2
Q = mp.mpf
rhs = (-21*L4(-2) - Q(65)/8*L4(Q(1)/3) + Q(13)/2*L4(Q(2)/3) + Q(137)/6*L4(Q(1)/2) + Q(39)/64*L4(Q(1)/9)
       + Q(15)/2*L3(Q(1)/3)*l2 + Q(37)/2*L3(-Q(1)/2)*l2 - Q(5)/4*L3(-Q(1)/2)*l3
       + Q(119)/16*z3*l2 - Q(37)/32*z3*l3
       + L2(Q(1)/3)**2 - Q(23)/4*L2(Q(1)/3)*l2**2 + Q(7)/4*L2(Q(1)/3)*l2*l3 + L2(Q(1)/3)*l3**2 - p2/6*L2(Q(1)/3)
       - Q(455)/144*p2*l2**2 + Q(9)/8*p2*l2*l3 - Q(7)/32*p2*l3**2 - Q(949)/2160*p2**2
       - Q(211)/144*l2**4 + Q(101)/24*l2**3*l3 - Q(15)/8*l2**2*l3**2 - Q(35)/24*l2*l3**3 + Q(113)/192*l3**4)
print('F4 =', mp.nstr(F4, 30)); print('F4 - closed form (rel1) =', mp.nstr(F4 - rhs, 5))
# The Li4(-2) etc. are real here? Li4(-2) real; fine.
