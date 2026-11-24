import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from flint import arb, acb, ctx
ctx.dps = 40
from mpmath import mp, mpf, pi, sqrt, exp, beta
mp.dps = 50
j = 1j
omega_m = exp(j*pi/4)
N = 400
y = sqrt(pi/(2*N)); s0 = omega_m/y

# reference h_j via quadrature (t-variable): h_j = omega^j H^{(j)}(s0)/(j! H(s0))
def Hc(s): return beta(s/2, mpf(3)/2)/2
def Hder(m):
    return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
H0 = Hc(s0)
h_ref = [omega_m**m*Hder(m)/H0/mp.factorial(m) for m in range(1,4)]
print("reference h1,h2,h3:", [mp.nstr(h, 10) for h in h_ref])

# polygamma version
PI = arb.pi(); SQ2 = arb(2).sqrt()
OMEGA = acb(1,1)/SQ2
s0a = OMEGA/arb(sqrt(float(pi/(2*N))))
Ls = []
for jj in range(1, 6):
    L1 = arb(2)**(-(jj)) * ((s0a/2).digamma() if jj==1 else (s0a/2).polygamma(jj-1))
    L2 = arb(2)**(-(jj)) * ((s0a/2+arb(3)/2).digamma() if jj==1 else (s0a/2+arb(3)/2).polygamma(jj-1))
    Ls.append(L1 - L2)
print("Ls:", [L.str(6) for L in Ls])
Ha = PI.sqrt()/2 * (s0a/2+1).gamma() * (s0a/2+arb(3)/2).rgamma()
Hd = [Ha]
for m in range(1, 4):
    Hd.append(sum(arb(arb(m-1).fac()/(arb(k).fac()*arb(m-1-k).fac()))*Ls[k]*Hd[m-1-k] for k in range(m)))
h1 = OMEGA*Hd[1]/Ha
h2 = OMEGA**2*Hd[2]/(2*Ha)
h3 = OMEGA**3*Hd[3]/(6*Ha)
print("polygamma h1,h2,h3:", [h.str(8) for h in [h1,h2,h3]])
