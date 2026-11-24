import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
from dcy2_common import df2, Hc, ell_coeff, omega, wb, j

N = mpf(400)
y = sqrt(pi/(2*N)); rho = 1/y
s0 = rho*omega
l1,l2,l3,l4,l5 = (ell_coeff(N,m) for m in range(1,6))
l6 = ell_coeff(N,6)
def Hder(m):
    return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
H0 = Hc(s0)
h1,h2,h3 = (Hder(m)/H0/mp.factorial(m) for m in range(1,4))
t0 = -l1/(2*l2)
def mu(p):
    base = exp(-l2*t0**2) * sqrt(-pi/l2)
    return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))

# exact ell(t)
def ell(t):
    return (N-1)*log(1+y*t) - N*log(1+y*t/(1+y*wb)) + j*pi*omega*t/2

# TRUE full-line integral of e^{ell} * F  (F = H-ratio)
I_FL = mp.quad(lambda t: exp(ell(t))*Hc(s0+omega*t)/H0, [-rho+mpf('0.001'), 0, 5, 20, mp.inf])

# moment sum with current P
P = {}
def add(p,c): P[p] = P.get(p,0)+c
add(0,1); add(1,h1); add(2,h2); add(3,h3)
add(3,l3); add(4,l4); add(5,l5)
add(6,l3**2/2); add(7,l3*l4); add(9,l3**3/6)
add(4,h1*l3); add(5,h1*l4); add(5,h2*l3); add(7,h1*l3**2/2)
M = sum(cp*mu(p) for p,cp in P.items())
print("I_FL (true, full integrand) =", I_FL)
print("M (moment sum)              =", M)
print("M - I_FL                    =", M - I_FL, " rel:", abs(M-I_FL)/abs(I_FL))
print("mu0 =", mu(0))
print("(M-I_FL)/y =", (M-I_FL)/y)
