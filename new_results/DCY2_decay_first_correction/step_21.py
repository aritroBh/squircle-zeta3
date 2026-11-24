from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
j = 1j
def Hc(s): return beta(s/2, mpf(3)/2)/2
N = mpf(10000)
y = sqrt(pi/(2*N)); v = sqrt(y); rho = 1/y
omega = exp(j*pi/4); wb = exp(-j*pi/4)
s0 = rho*omega
def phi(s): return (N-1)*log(s) - N*log(1+s) + j*pi*s/2

b = pi*wb/2
def G1(tau): return -(pi*wb**2/2)*tau + (pi*wb/2)*tau**3
def G2(tau): return (3*pi*wb**2/4)*tau**2 - (pi*wb/2)*tau**4
def Gser(tau): return -b*tau**2 + v*G1(tau) + v**2*G2(tau)
for tau in [mpf('0.3'), mpf('0.7'), mpf('1.5')]:
    t = tau/sqrt(y)
    Ge = phi(s0+omega*t)-phi(s0)
    print(f"tau={tau}: exact={Ge}, ser={Gser(tau)}, diff={Ge-Gser(tau)}")
    Hr = Hc(s0+omega*t)/Hc(s0)
    Hser = 1 - mpf(3)/2*v*tau + mpf(15)/8*v**2*tau**2
    print(f"   Hratio exact={Hr}, ser={Hser}, diff={Hr-Hser}")
