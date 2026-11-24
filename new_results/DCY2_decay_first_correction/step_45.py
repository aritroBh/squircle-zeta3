import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp, gamma, log, factorial
mp.dps = 50
j = 1j
omega = exp(j*pi/4); wb = exp(-j*pi/4)
bR = pi/(2*sqrt(2))
N = mpf(10000)
y = sqrt(pi/(2*N)); rho = 1/y
kap0 = 2.48
kap = sqrt(kap0**2 + 4*log(N/400))
def ell(m):
    if m == 1: return (N-1)*y - N*y/(1+y*wb) + j*pi*omega/2
    return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
l = {m: ell(m) for m in range(1,7)}
c2m = 0.9*bR*y
sig = 1/sqrt(2*c2m)
Delta = kap*sig
yD = y*Delta
L3,L4,L5 = abs(l[3]),abs(l[4]),abs(l[5])
B3 = (N*y**4 + y**3/3)/(1-yD); B3sup = B3*Delta**3
B6 = (N*y**7 + y**6/6)/(1-yD); B6sup = B6*Delta**6
alpha = L3 + L4*Delta + L5*Delta**2
beta2 = L3*L5 + L4**2/2 + L4*L5*Delta + L5**2*Delta**2/2
beta3 = (3*L3**2*L4 + (3*L3**2*L5+3*L3*L4**2)*Delta + (6*L3*L4*L5+L4**3)*Delta**2
         + (3*L3*L5**2+3*L4**2*L5)*Delta**3 + 3*L4*L5**2*Delta**4 + L5**3*Delta**5)/6
print(f"kap={float(kap):.3f} Delta={float(Delta):.2f} yD={float(yD):.3f}")
print(f"B3sup={float(B3sup):.3f} eB3={float(exp(B3sup)):.2f}  alpha*Delta^3={float(alpha*Delta**3):.3f} eA={float(exp(alpha*Delta**3)):.2f}")
x = 1/(sqrt(2)*y)
M4num = sqrt(2)*gamma(mpf(11)/2)*x**(-mpf(11)/2)
Cg=mpf(1); g_lo = 1-3*y/4-Cg*y**2
H_lo = sqrt(pi/2)*rho**(-mpf(3)/2)*g_lo
M4 = M4num/H_lo
ReL1 = abs(l[1].real); tR = ReL1/(2*c2m); efac = exp(ReL1**2/(4*c2m))
def Mtilde(jj):
    return efac*sum(mp.binomial(jj,i)*tR**(jj-i)*c2m**(-(i+1)/2)*gamma((i+1)/2) for i in range(jj+1))
for jj in [4,6,8,10,12]:
    print(f"Mtilde({jj}) = {float(Mtilde(jj)):.3e}")
print(f"M4/24*Mtilde(4) = {float(M4/24*Mtilde(4)):.3e}")
print(f"beta2*Mtilde(8) = {float(beta2*Mtilde(8)):.3e}")
print(f"beta3*Mtilde(10) = {float(beta3*Mtilde(10)):.3e}")
print(f"alpha^4/24*eA*Mtilde(12) = {float(alpha**4/24*exp(alpha*Delta**3)*Mtilde(12)):.3e}")
print(f"B6*eB6*Mtilde(6) = {float(B6*exp(B6sup)*Mtilde(6)):.3e}")
