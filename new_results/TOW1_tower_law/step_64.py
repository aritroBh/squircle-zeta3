from towlib import *
import sympy as sp
# numeric values at p=3, key (3,2), d=5
r0 = rho_c((3,2),0)
def sy(q): return q.sympy()
print('rho(2)/rho:', sy(rho_c((3,2),2)/r0))
print('R2(5,3) =', sp.Rational(1,2)*( (15-3)/2 - sp.I*(2*(5-6)-10*sp.I-4*4*(5-6)/7)/2 - 9 ))
print('rB0 exact rho_43(0)/rho_32(0):', sy(rho_c((4,3),0)/r0))
print('rB2 exact rho_43(2)/rho_32(0):', sy(rho_c((4,3),2)/r0))
# ehat values exact
print('ehat33(0,0)/rho:', sy(ehat_pp(3,0,0)/r0))
print('ehat33(1,0)/rho:', sy(ehat_pp(3,1,0)/r0))
print('ehat33(2,0)/rho:', sy(ehat_pp(3,2,0)/r0))
print('ehat33(2,1)/rho:', sy(ehat_pp(3,2,1)/r0))
# exact feeds at delta=2 for key (3,2) from towlib recursion pieces:
# DL1 from (4,2): (5)*sum tau1_{2-dp} e_42(dp)/(10-dp-1)... compute total RHS and LHS
lhs = rho_c((3,2),2) + e_off((3,2),2,1)
rhs = LambdaE_ij(3,2,2,0)
print('exact lhs==rhs:', lhs==rhs, sy((lhs-rhs)/r0))
# my fold formula pieces at p=3: m=6, g_6,g_7,g_8
ed = {0: ehat_pp(3,0,0), 1: ehat_pp(3,1,0), 2: ehat_pp(3,2,0)-ehat_pp(3,2,1)}
for mu in (6,7,8):
    g = ZERO
    for dp in range(0, 4*3-4-mu+1):
        g = g - ed[dp] * tau1(4*3-4-mu-dp).scale(Fr(1, 4*3-3-dp))
    print(f'mu={mu}: g/rho =', sy(g/r0))
