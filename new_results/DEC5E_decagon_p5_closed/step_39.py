import mpmath as mp
mp.mp.dps = 800
phi = (1+mp.sqrt(5))/2
lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi**2
# Lewin identity (from Lewin 3-term Landen at x=1/phi^2 + duplication):
L34 = mp.polylog(3, 1/phi**2)
rhs = mp.mpf(4)/5*z3 - mp.mpf(2)/15*pi2*lf + mp.mpf(2)/3*lf**3
print("Li3(1/phi^2) = 4/5 z3 - 2/15 pi^2 lf + 2/3 lf^3:", mp.nstr(abs(L34-rhs), 3))
r2 = mp.polylog(3,1/phi) + mp.polylog(3,-1/phi) - (z3/5 - pi2*lf/30 + lf**3/6)
print("Li3(1/phi)+Li3(-1/phi) = 1/5 z3 - 1/30 pi^2 lf + 1/6 lf^3:", mp.nstr(abs(r2), 3))
# joint relations among Li4 points + ring?
ring4 = [mp.pi**4, z3*lf, z3*mp.log(2), pi2*lf**2, pi2*lf*mp.log(2), pi2*mp.log(2)**2, lf**4, lf**3*mp.log(2), lf**2*mp.log(2)**2, lf*mp.log(2)**3, mp.log(2)**4]
rel = mp.pslq(ring4 + [mp.polylog(4,-1/phi), mp.polylog(4,1/phi**2)], tol=1e-700, maxcoeff=10**15, maxsteps=10000)
print("Li4 pair + ring4 relation:", rel)
ring5 = [z5, pi2*z3, z3*lf**2, z3*lf*mp.log(2), z3*mp.log(2)**2, mp.pi**4*lf, mp.pi**4*mp.log(2)]
rel5 = mp.pslq(ring5 + [mp.polylog(5,-1/phi), mp.polylog(5,1/phi**2)], tol=1e-700, maxcoeff=10**15, maxsteps=10000)
print("Li5 pair + ring5(partial) relation:", rel5)
