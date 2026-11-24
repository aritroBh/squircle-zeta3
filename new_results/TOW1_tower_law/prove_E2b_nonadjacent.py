from proof_common import *
d, p = sp.symbols('d p', integer=True, positive=True)
# E2b (depth-2 value identity), non-adjacent key (p,q), q=d-p, normalized by rho_{p,q}(0).
LHS = R2(d,p) - 1/(2*d-3)
rA0 = i*(2*d-1)/(p+1)
rB0 = -i*(2*d-1)/(d-p+1)
RHS = ( (p+1)*( t[(1,2)]*rA0/(2*d-1) + t[(1,1)]*(-i*(d-1))*rA0/(2*d-2) + t[(1,0)]*(R2(d+1,p+1)*rA0)/(2*d-3) )
      + (d-p+1)*( t[(2,2)]*rB0/(2*d-1) + t[(2,1)]*(-i*(d-1))*rB0/(2*d-2) + t[(2,0)]*(R2(d+1,p)*rB0)/(2*d-3) ) )
print('E2b-non-adjacent remainder:', sp.cancel(sp.together(LHS - RHS)))
