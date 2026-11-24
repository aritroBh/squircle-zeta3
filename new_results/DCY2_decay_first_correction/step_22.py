from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
j = 1j
def Hc(s): return beta(s/2, mpf(3)/2)/2
for N in [mpf(10000), mpf(1000000)]:
    y = sqrt(pi/(2*N)); v = sqrt(y); rho = 1/y
    omega = exp(j*pi/4); wb = exp(-j*pi/4)
    s0 = rho*omega
    def phi(s): return (N-1)*log(s) - N*log(1+s) + j*pi*s/2
    b = pi*wb/2
    def G1(tau): return -(pi*wb**2/2)*tau + (pi*wb/2)*tau**3
    tau = mpf('0.7')
    t = tau/sqrt(y)
    Ge = phi(s0+omega*t)-phi(s0)
    G2_fit = (Ge + b*tau**2 - v*G1(tau))/v**2
    print(f"N={N}: G2_fit = {G2_fit}")
print("my G2(0.7) =", (3*pi*exp(-j*pi/2)/4)*mpf('0.49') - (pi*exp(-j*pi/4)/2)*mpf('0.7')**4)
# Also check G1
N = mpf(1000000)
y = sqrt(pi/(2*N)); v = sqrt(y); rho=1/y; omega=exp(j*pi/4); wb=exp(-j*pi/4); s0=rho*omega
def phi(s): return (N-1)*log(s) - N*log(1+s) + j*pi*s/2
b = pi*wb/2
tau = mpf('0.7'); t = tau/sqrt(y)
Ge = phi(s0+omega*t)-phi(s0)
print("G1_fit =", (Ge + b*tau**2)/v)
print("my G1(0.7) =", -(pi*wb**2/2)*tau + (pi*wb/2)*tau**3)
