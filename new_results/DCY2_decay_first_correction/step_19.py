from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
j = 1j

def Hc(s): return beta(s/2, mpf(3)/2)/2

for N in [mpf(400), mpf(10000)]:
    y = sqrt(pi/(2*N)); rho = 1/y
    omega = exp(j*pi/4); wb = exp(-j*pi/4)
    s0 = rho*omega
    theta = sqrt(pi*N)
    def phi(s): return (N-1)*log(s) - N*log(1+s) + j*pi*s/2
    tmax = 1/sqrt(y)  # tau in [-1/sqrt(y), inf)
    def integr(tau):
        t = tau/sqrt(y)
        return exp(phi(s0+omega*t)-phi(s0)) * Hc(s0+omega*t)/Hc(s0)
    Jnum = mp.quad(integr, [-tmax+mpf('0.001'), 0, 5, mp.inf])
    b = pi*wb/2
    c1 = (1+1j)*(15-pi**2)/(8*sqrt(2)*pi)
    Jformal = sqrt(pi/b)*(1 + c1*y)
    print(f"N={N}: Jnum = {Jnum}")
    print(f"   Jformal = {Jformal}, ratio = {Jnum/Jformal}")
