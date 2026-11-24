import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
j = 1j
omega = exp(j*pi/4)
def Hc(s): return beta(s/2, mpf(3)/2)/2
def df2(k):
    r = mpf(1)
    while k > 1: r *= k; k -= 2
    return r
varpi1 = ((24+12*pi+pi**2) + 1j*(24-12*pi+pi**2))/(48*sqrt(pi))
def model(N):
    y = sqrt(pi/(2*N)); rho = 1/y
    wb = exp(-j*pi/4); s0 = rho*omega
    phi_s0 = (N-1)*log(s0) - N*log(1+s0) + j*pi*s0/2
    def ell(m):
        if m == 1: return (N-1)*y - N*y/(1+y*wb) + j*pi*omega/2
        return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
    l1,l2,l3,l4,l5 = (ell(m) for m in range(1,6))
    def Hder(m):
        return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
    H0 = Hc(s0)
    h1,h2,h3 = (omega**m*Hder(m)/H0/mp.factorial(m) for m in range(1,4))
    t0 = -l1/(2*l2)
    def mu(p):
        base = exp(-l2*t0**2) * sqrt(-pi/l2)
        return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))
    P = {}
    def add(p,c): P[p] = P.get(p,0)+c
    add(0,1); add(1,h1); add(2,h2); add(3,h3)
    add(3,l3); add(4,l4); add(5,l5)
    add(6,l3**2/2); add(7,l3*l4); add(9,l3**3/6)
    add(4,h1*l3); add(5,h1*l4); add(5,h2*l3); add(7,h1*l3**2/2)
    return (2/pi)*exp(phi_s0)*H0*omega*sum(cp*mu(p) for p,cp in P.items())

def target(N):
    theta = sqrt(pi*N)
    return -1j*sqrt(pi)/N*exp(-theta+j*theta)*(1+varpi1/sqrt(N))
def leading(N):
    theta = sqrt(pi*N)
    return sqrt(pi)/N*exp(-theta)

for N in [400, 2500, 10000]:
    N = mpf(N)
    Mo = model(N); T1 = target(N); L = leading(N)
    print(f"N={int(N)}  |model-T1|/(L/N) = {float(abs(Mo-T1)/(L/N)):.5f}")
