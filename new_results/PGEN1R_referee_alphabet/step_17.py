import mpmath as mp
mp.mp.dps = 400
phi = (1+mp.sqrt(5))/2
s5 = mp.sqrt(5)
z3 = mp.zeta(3)
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
# L(chi5,3) from Hurwitz zeta
L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
# cross-check with direct sum
L5b = mp.nsum(lambda k: mp.kronecker(5,k)/mp.mpf(k)**3, [1, mp.inf], tol=mp.mpf('1e-380'))
print('L5 =', mp.nstr(L5, 50))
print('L5 direct diff:', mp.nstr(L5-L5b, 5))
# golden chi3 identity D
chi3 = (Li3(1/phi)-Li3(-1/phi))/2
D_id = chi3 - ((25*s5/48)*L5 - (mp.pi**2/12)*mp.log(phi) + mp.mpf(1)/4*mp.log(phi)**3)
print('claim D residual (400 dps):', mp.nstr(D_id, 5))
# controls: Li2 golden identities
c_Li2a = Li2(1/phi) - (mp.pi**2/10 - mp.log(phi)**2)
c_Li2b = Li2(-1/phi) - (-mp.pi**2/15 + mp.log(phi)**2/2)
print('Li2(1/phi) id residual:', mp.nstr(c_Li2a, 5))
print('Li2(-1/phi) id residual:', mp.nstr(c_Li2b, 5))
# Lewin Li3(1/phi^2): identify coefficients of [z3, pi^2 log phi, log^3 phi]
t = Li3(1/phi**2)
rel = mp.pslq([t, z3, mp.pi**2*mp.log(phi), mp.log(phi)**3], tol=mp.mpf('1e-350'), maxcoeff=10**6)
print('Lewin Li3(1/phi^2) pslq [t, z3, pi2logphi, log3phi]:', rel)
