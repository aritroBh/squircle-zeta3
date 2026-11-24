import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
from dcy2_common import df2, Hc, ell_coeff, omega, wb, j
N = mpf(400); y = sqrt(pi/(2*N)); s0 = omega/y
l1,l2,l3 = (ell_coeff(N,m) for m in range(1,4))
H0 = Hc(s0)
def ell(t):
    return (N-1)*log(1+y*t) - N*log(1+y*t/(1+y*wb)) + j*pi*omega*t/2
def R3(t): return ell(t) - l1*t - l2*t**2
def Ffull(t): return Hc(s0+omega*t)/H0 * exp(R3(t))
for t in [mpf('0.001'), mpf('-0.001'), mpf('0.01'), mpf('0.05')]:
    print(t, Ffull(t))
print("finite diff t1:", (Ffull(mpf('0.001'))-Ffull(-mpf('0.001')))/(2*mpf('0.001')))
print("h1 = -0.0663096769402 + 0.0636951459622j")
# also check R3 expansion: R3(t) should be ~ l3 t^3
print("R3(0.01) =", R3(mpf('0.01')), " l3*0.01^3 =", l3*mpf('0.01')**3)
