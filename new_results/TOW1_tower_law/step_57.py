from towlib import *
# exact pieces at p=2, key (2,1), d=3, normalize by rho_21(0)
r0 = rho_c((2,1),0)
print('rho_21(0) =', r0)
print('ehat00 (code) =', ehat_pp(2,0,0)/r0, ' formula t10*(4p-3)/p =', QSI(0,0,Fr(-1,2),Fr(-1,2)) if False else None)
import sympy as sp
t10v = (-1-sp.I)/2
print('formula ehat00/rho =', t10v*(4*2-3)/2)  # expect (ehat_22(0,0))/rho_21(0)
print('exact ehat22(0,0)/rho21(0) =', (ehat_pp(2,0,0)/r0).sympy())
print('exact ehat22(1,0)/rho21(0) =', (ehat_pp(2,1,0)/r0).sympy())
# formula rB0, rB1
print('rB0 formula =', (2*3-1)*(2*3+1)/((2+1)*2), ' exact rho_32(0)/rho_21(0) =', (rho_c((3,2),0)/r0).sympy())
print('rB1 exact =', (rho_c((3,2),1)/r0).sympy(), ' formula -i*d*rB0 =', -sp.I*3*sp.Rational(35,6))
# the fold feed at delta=1 from towlib directly
i0,j0 = 2,1
fold_part = Lambda_ij(i0,j0,4,1, lambda kk,dp: ZERO)  # only fold survives? no: Lambda has all feeds zeroed
print('fold exact (from LambdaE with all parents zero):', fold_part)
# better: compute fold via the e-families
def Udiag(dp):
    if dp==0: return ehat_pp(2,0,0)
    if dp==1: return ehat_pp(2,1,0)
    return ZERO
from math import comb as _c
m = 4*2-4-1
ft = ZERO
for mu in range(m, 4*2-3):
    g = ZERO
    for dp in range(0, 4*2-4-mu+1):
        g = g - Udiag(dp) * tau1(4*2-4-mu-dp).scale(Fr(1, 4*2-3-dp))
    ft = ft + g.conj_i() * (IQ**(-mu)) * Fr(_c(mu-2, mu-m)) * (IQ**(m-2))
ft = ft.scale(-2)
print('fold exact/rho =', (ft/r0).sympy())
print('feed_DL2 exact/rho =', ((2*( tau2(0)*ehat_pp(2,1,0).scale(Fr(1,2*3-2)) + tau2(1)*ehat_pp(2,0,0).scale(Fr(1,2*3-1)) ))/r0).sympy())
