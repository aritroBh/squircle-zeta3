import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
from dcy2_common import df2, Hc, ell_coeff, omega, wb, j

N = mpf(400)
y = sqrt(pi/(2*N)); rho = 1/y
s0 = rho*omega
l1,l2,l3,l4,l5,l6 = (ell_coeff(N,m) for m in range(1,7))
def Hder(m):
    return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
H0 = Hc(s0)
hs = [mpf(1)] + [Hder(m)/H0/mp.factorial(m) for m in range(1,7)]
def ell(t):
    return (N-1)*log(1+y*t) - N*log(1+y*t/(1+y*wb)) + j*pi*omega*t/2
def R3(t): return ell(t) - l1*t - l2*t**2
def Ffull(t): return Hc(s0+omega*t)/H0 * exp(R3(t))

# Taylor coefficients of Ffull at 0 by fitting values at small t (Vandermonde, t_k = k*eps)
eps = mpf('0.01')
ns = 8
A = []
Bv = []
for k in range(1, ns+1):
    t = eps*k
    A.append([t**p for p in range(ns)])
    Bv.append(Ffull(t))
coef = mp.lu_solve(mp.matrix(A), mp.matrix(Bv))
print("fitted Taylor coeffs of Ffull:")
for p in range(1,7):
    print(f"  t^{p}: {coef[p]}")
print()
print("h coefficients:", [mp.nstr(h,10) for h in hs])
print("P expects: t^3:", hs[3]+l3, " t^4:", l4+hs[1]*l3, " t^5:", l5+hs[1]*l4+hs[2]*l3, " t^6:", l3**2/2)
