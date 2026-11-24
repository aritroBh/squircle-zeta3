from towlib import *
from math import comb as _c
p0 = 2
def Udiag(dp):
    return ehat_pp(p0,dp,0)
for mu in (3,4):
    g = ZERO
    for dp in range(0, 4*p0-4-mu+1):
        g = g - Udiag(dp) * tau1(4*p0-4-mu-dp).scale(Fr(1, 4*p0-3-dp))
    print(f'mu={mu}: g = {g}, conj = {g.conj_i()}')
# the full fold from Lambda_ij (single out the fold by zeroing other parents)
fold_only = Lambda_ij(2,1,4,1, lambda kk,dp: Udiag(dp) if kk==(2,2) else ZERO)
print('fold via Lambda =', fold_only)
# subtract the DL2-direct part (which also uses (2,2))
dl2 = ( tau2(0)*ehat_pp(2,1,0).scale(Fr(2, 2*3-2)) + tau2(1)*ehat_pp(2,0,0).scale(Fr(2, 2*3-1)) )
print('DL2-direct part =', dl2)
print('fold proper =', fold_only - dl2)
print('fold proper/rho =', ((fold_only-dl2)/rho_c((2,1),0)).sympy())
