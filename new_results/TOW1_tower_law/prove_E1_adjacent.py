from proof_common import *
p = sp.Symbol('p', integer=True, positive=True)
d = 2*p-1
# E1 (depth-1 identity), adjacent key (p,p-1), normalized by rho_{p,p-1}(0).
LHS = -i*(d-2)
rA0 = i*(2*d-1)/(p+1)                       # rho_{p+1,p-1}(0)/rho_{p,p-1}(0)
rA1 = -i*(d-1)*rA0
rB0 = (2*d-1)*(2*d+1)/((p+1)*p)             # rho_{p+1,p}(0)/rho_{p,p-1}(0)
rB1 = -i*d*rB0
ehat00 = t[(1,0)]*(4*p-3)/p                 # ehat_pp(0,0)/rho = h_hat_p/rho
ehat10 = (p+1)*( t[(1,0)]*rB1/(4*p-2) + t[(1,1)]*rB0/(4*p-1) )
feed_DL1 = (p+1)*( t[(1,0)]*rA1/(2*d-2) + t[(1,1)]*rA0/(2*d-1) )
feed_DL2 = p*( t[(2,0)]*ehat10/(2*d-2) + t[(2,1)]*ehat00/(2*d-1) )
ehat00c = sp.conjugate(ehat00)
ehat10c = (p+1)*( t[(2,0)]*sp.conjugate(rB1)/(4*p-2) + (-t[(2,1)])*sp.conjugate(rB0)/(4*p-1) )
g_m1c = -( ehat00c*(-t[(2,1)])/(4*p-3) + ehat10c*t[(2,0)]/(4*p-4) )
g_m0c = -( ehat00c*t[(2,0)]/(4*p-3) )
feed_fold = p*( g_m1c - i*(4*p-6)*g_m0c )
RHS = feed_DL1 + feed_DL2 + feed_fold
print('E1-adjacent remainder:', sp.cancel(sp.together(LHS - RHS)))
