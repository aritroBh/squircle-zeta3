from towlib import *
print('rho_32:', [rho_c((3,2),dp) for dp in (0,1,2)])
print('u_32(5):', [u(3,2,5,dp) for dp in (0,1,2)])
lhs = u(2,2,5,2)
rhs = ( e_off((3,2),2,0) * tau1(0).scale(Fr(3, 5))
      + e_off((3,2),1,0) * tau1(1).scale(Fr(3, 6))
      + e_off((3,2),0,0) * tau1(2).scale(Fr(3, 7)) )
print('u_22(5,2) =', lhs)
print('feed      =', rhs)
print('u_32(5,0..2) == rho_32:', [u(3,2,5,dp)==rho_c((3,2),dp) for dp in (0,1,2)])
# direct Lambda check
print('Lambda:', Lambda_ij(2,2,5,2, lambda kk,dp: u(kk[0],kk[1],5,dp)))
