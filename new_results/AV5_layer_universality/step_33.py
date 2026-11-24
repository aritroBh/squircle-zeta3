import mpmath
from mpmath import mp, mpf, zeta, log, nstr
mp.dps = 90
L = log(2); Z2 = zeta(2); Z3 = zeta(3)
# T = Z3/8, basis {L^3, Z2 L, Z3}
T = Z3/8
basis = [L**3, Z2*L, Z3, T]
rel = mpmath.pslq(basis)
print("default pslq:", rel)
rel2 = mpmath.pslq(basis, tol=mpf(10)**(-70), maxcoeff=10**6, maxsteps=10000)
print("strict pslq:", rel2)
