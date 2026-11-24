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
h1,h2,h3 = (Hder(m)/H0/mp.factorial(m) for m in range(1,4))
def ell(t):
    return (N-1)*log(1+y*t) - N*log(1+y*t/(1+y*wb)) + j*pi*omega*t/2
def R3(t): return ell(t) - l1*t - l2*t**2
def Ffull(t): return Hc(s0+omega*t)/H0 * exp(R3(t))

P = {}
def add(p,c): P[p] = P.get(p,0)+c
add(0,1); add(1,h1); add(2,h2); add(3,h3)
add(3,l3); add(4,l4); add(5,l5)
add(6,l3**2/2); add(7,l3*l4); add(9,l3**3/6)
add(4,h1*l3); add(5,h1*l4); add(5,h2*l3); add(7,h1*l3**2/2)
def Pval(t): return sum(cp*t**p for p,cp in P.items())

for t in [mpf('0.5'), mpf('1.0'), mpf('2.0')]:
    print(f"t={t}: Ffull={Ffull(t)}")
    print(f"      P    ={Pval(t)}   diff={Ffull(t)-Pval(t)}  (diff)/y = {(Ffull(t)-Pval(t))/y}")
