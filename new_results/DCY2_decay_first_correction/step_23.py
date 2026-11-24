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
    tau = mpf('0.7'); t = tau/sqrt(y)
    Ge = phi(s0+omega*t)-phi(s0)
    G1_fit = (Ge + b*tau**2)/v
    print(f"N={N}: G1_fit={G1_fit}")
# my corrected G1: cubic (pi wb/2) tau^3, linear: -(pi wb^2/2) tau - tau
wb = exp(-j*pi/4)
tau = mpf('0.7')
print("corrected G1 =", (pi*wb/2)*tau**3 + (-(pi*wb**2/2)-1)*tau)
# G2 candidates: quadratic: 3 pi wb^2 tau^2 /4 + tau^2/2 ; quartic: -(pi wb/2) tau^4
print("corrected G2 =", (3*pi*wb**2/4+mpf(1)/2)*tau**2 - (pi*wb/2)*tau**4)
