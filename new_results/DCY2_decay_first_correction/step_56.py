import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
j = 1j
omega = exp(j*pi/4); wb = exp(-j*pi/4)
b = pi*wb/2
def Hc(s): return beta(s/2, mpf(3)/2)/2
def df2(k):
    r = mpf(1)
    while k > 1: r *= k; k -= 2
    return r
N = mpf(10**7)
y = sqrt(pi/(2*N)); s0 = omega/y
def ell(m):
    if m == 1: return (N-1)*y - N*y/(1+y*wb) + j*pi*omega/2
    return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
l1,l2,l3,l4,l5 = (ell(m) for m in range(1,6))
t0 = -l1/(2*l2)
def mu(p):
    base = exp(-l2*t0**2)*sqrt(-pi/l2)
    return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))
def Hder(m):
    return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
H0 = Hc(s0)
h1,h2,h3 = (omega**m*Hder(m)/H0/mp.factorial(m) for m in range(1,4))
mu0 = mu(0)
print("l4*M4/y    =", l4*mu(4)/mu0/y, " expect 3om/2pi =", 3*omega/(2*pi))
print("h1l3*M4/y  =", h1*l3*mu(4)/mu0/y, " expect -9om/4pi =", -9*omega/(4*pi))
# check M4 itself
print("M4 =", mu(4)/mu0, " 3/(2by)^2 =", 3/(2*b*y)**2, " 6t0^2/(2by) =", 6*t0**2/(2*b*y))
# l4 expansion check
print("l4 vs pi wb y^3/2:", l4, pi*wb*y**3/2)
print("h1 vs -3y/2:", h1, -3*y/2)
print("l3 vs pi wb y^2/2:", l3, pi*wb*y**2/2)
